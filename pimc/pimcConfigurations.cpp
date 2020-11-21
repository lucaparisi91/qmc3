#include "pimcConfigurations.h"
#include "qmcExceptions.h"
#include <filesystem>
#include <sstream>
#include <algorithm>


namespace fs = std::filesystem;

namespace pimc
{

int getTotalSize( const std::vector<particleGroup> & groups )
{
    int nMax=0;
    for (const auto & group : groups)
    {
        nMax=std::max(group.iEndExtended,nMax);
    }

    return nMax+1;
}

int getNParticles( const std::vector<particleGroup> & groups )
{
    int N=0;
    for (const auto & group : groups)
    {
        N+=group.iEnd - group.iStart + 1;
    }

    return N;
}

chain::chain() :
prev(-1) , next(-1) , head(0) , tail(-1)
    {

    }


void pimcConfigurations::setHead( int iChain, int newHead )
{
    int oldHead=_chains[iChain].head;

    int delta = newHead - oldHead;

    if ( delta > 0)
    {
        createBeads({oldHead,newHead-1},iChain);
    }
    else
    {
        deleteBeads({newHead,oldHead},iChain);
    }
    
    auto & currentChain=_chains[iChain];
    int oldNext=currentChain.next;

    currentChain.head=newHead;
    currentChain.next=-1;

    if ( oldNext  >= 0 )
    {
        setTail(oldNext,_chains[oldNext].tail);
        _heads.push_back( iChain);
    }

    
};


void pimcConfigurations::setTail( int iChain, int newTail )
{
    int oldTail=_chains[iChain].tail;

    int delta = newTail - oldTail;

    if ( delta < 0)
    {
        createBeads({newTail+1,oldTail},iChain);
    }
    else
    {
        deleteBeads({oldTail+1,newTail},iChain);
    }

    auto & currentChain=_chains[iChain];

    int oldPrev=currentChain.prev;

    currentChain.tail=newTail;
    currentChain.prev=-1;


    if ( oldPrev  >= 0 )
    {
        setHead(oldPrev,_chains[oldPrev].head);
        _tails.push_back( iChain);
    }

};


pimcConfigurations::pimcConfigurations(
    size_t timeSlices, int dimensions, const std::vector<particleGroup> &  particleGroups_) :
    particleGroups(particleGroups_),
    M(timeSlices),
    N(getTotalSize(particleGroups_)), // number of chains(includign padding chains)
    _data(N,dimensions,2*timeSlices ),// contains a copy for buffer operations
     _mask( timeSlices+1,N),
    _nParticles(getNParticles(particleGroups_))// number of chains without the padding
{
    _chains.resize(N);
    for(int i=0;i<N;i++)
    {
        _chains[i].tail=-1;
        _chains[i].head=nBeads();
        _chains[i].next=i;
        _chains[i].prev=i;
    }

}; 

void pimcConfigurations::fillHead(int iChain)
{
    const auto & currentChain = _chains[iChain];

    if ( currentChain.next != -1 )
    {

        for (int d=0;d<getDimensions();d++)
        {
            _data(iChain,d,currentChain.head )=_data(currentChain.next,d, 
            _chains[currentChain.next].tail + 1
            );
        }

    }

}

void pimcConfigurations::fillHeads()
{
   for( const auto & group : particleGroups )
   {
       for(int i=group.iStart;i<=group.iEnd;i++)
       {
            fillHead(i);
       }
       
   }    

}



void pimcConfigurations::deleteBeads(   std::array<int,2> timeRange, int iChain )
{
    for (int t=timeRange[0] ; t<= timeRange[1] ;t++   )
    {
        _mask(t,iChain)=0;
    }
}

void pimcConfigurations::createBeads(   std::array<int,2> timeRange, int iChain )
{
    for (int t=timeRange[0] ; t<= timeRange[1] ;t++   )
    {
        _mask(t,iChain)=1;
    }
}

int pimcConfigurations::pushChain( particleGroup & group)
{
    group.iEnd+=1;

    if (group.iEnd > group.iEndExtended )
    {
        throw invalidState("No availible memory left for this group. ");
    }

    if ( _chains[group.iEnd].hasHead() )
    {
        _heads.push_back(group.iEnd);
    }

    if ( _chains[group.iEnd].hasTail() )
    {
        _tails.push_back(group.iEnd);
    }

    return group.iEnd;
}


int pimcConfigurations::pushChain( int iGroup)
{
    auto & group = particleGroups[iGroup];
    return pushChain(group);
}


void pimcConfigurations::removeChain( int iChain)
{
    auto & group=getGroupByChain(iChain);

    // unregister chain
    if ( !(getChain(iChain).hasHead()) or !( getChain(iChain).hasTail()) )
    {
        throw invalidState("Cannot remove chain with no head or tail");
    }

    swap(iChain,group.iEnd);

    deleteTailFromList(group.iEnd);
    deleteHeadFromList(group.iEnd);


    group.iEnd-=1;

}

void pimcConfigurations::deleteHeadFromList(int iChain)
{
    // delete the head from the list of heads
    auto it = std::find(_heads.begin(),_heads.end(),iChain) ;
    if ( it== _heads.end() )
    {
        throw invalidState("Could not find the head among registered heads.");
    }

    std::swap(*it,*(_heads.end() - 1) );
    _heads.resize(_heads.size() - 1);

}

void pimcConfigurations::deleteTailFromList(int iChain)
{
    // delete the tail from the list of tails
    auto it = std::find(_tails.begin(),_tails.end(),iChain) ;
    if ( it== _tails.end() )
    {
        throw invalidState("Could not find the tail among registered tails.");
    }
    std::swap(*it,*(_tails.end() - 1) );
    _tails.resize(_tails.size() - 1);


}



void pimcConfigurations::join( int iChainLeft, int iChainRight)
{
    setHead(iChainLeft,_chains[iChainLeft].head);
    setTail(iChainRight,_chains[iChainRight].tail);

    

    
    _chains[iChainLeft].next=iChainRight;
    _chains[iChainRight].prev=iChainLeft;

    deleteHeadFromList(iChainLeft);
    deleteTailFromList(iChainRight);
    
    

}



void pimcConfigurations::swapTails(int iChain1, int iChain2)
{
    auto headChain=getChain(iChain1);
    auto partnerChain = getChain(iChain2);

    setTail(iChain2,headChain.tail);
    setTail(iChain1,partnerChain.tail);

        if ( partnerChain.prev != -1)
        {
            join(partnerChain.prev,headChain.tail);
        }
        if ( headChain.prev != -1)
        {
            join(headChain.prev,headChain.next);
        }
}

void pimcConfigurations::save(const std::string & dirname,const std::string & format ) const
{

    if ( ! fs::exists(dirname) ) 
    { 
        fs::create_directory(dirname); // create src folder
    }

    std::ofstream f;

    if ( format == "csv" )
    {
    
        f.open(dirname + "/particles.dat");
        //f.write(reinterpret_cast<char*>(_data.data()), nBeads()*nChains()*getDimensions()*sizeof(double));
        f << "particle time x y z" << std::endl;
        f << std::setprecision(7);

        for (int t=0;t<nBeads() ;t++ )
        {
                for (int i=0;i<nChains();i++)
                {
                    if (_mask(t,i) == 1)
                        {
                            f<< i << " " << t << " ";
                            for (int d=0;d<getDimensions();d++)
                            {
                                
                                f <<  _data(i,d,t) << " ";
                            }
                            f<< std::endl;
                        }
                    
                }
        }
        f.close();
    }
    else if (format == "binary")
    {
        f.open(dirname + "/particles.dat",std::ios::binary);
        f.write(reinterpret_cast<const char*>(_data.data()), nBeads()*nChains()*getDimensions()*sizeof(double));
    }


    nlohmann::json j;

    j["nChains"]= nChains();
    j["dimensions"]=getDimensions();
    j["timeSlices"]=nBeads();
    j["format"]=format;

    std::vector<nlohmann::json> jGroups;
    for( const auto & group  : particleGroups)
    {
        nlohmann::json jGroup;
        jGroup["iStart"]=group.iStart;
        jGroup["iEnd"]=group.iEnd;
        jGroup["iEndExtended"]=group.iEndExtended;
        jGroups.push_back(jGroup);
    }

    j["groups"]=jGroups;
    f.open(dirname + "/description.dat",std::ios::binary);
    f << j ;
    f.close();

    throw missingImplementation("Saving of chain info not yet implemented");

}

void pimcConfigurations::load(const std::string & dirname)
{

    if ( ! fs::exists(dirname) ) 
    { 
        throw invalidInput("Folder " + dirname + " does not exist.");
    }

    std::ifstream f;
    f.open(dirname + "/description.dat");
    nlohmann::json j;
    
    f >> j;

    f.close();


    N = j["nChains"].get<int>();
    M=j["timeSlices"].get<int>();

    _data.resize(N,getDimensions(),M);
    _mask=mask(nBeads(),nChains()) ;
    _mask.setConstant(0);

    auto format = j["format"].get<std::string>();

    if (format == "csv")
    {
        std::string dummy;
        const int nFields = 5;
        f.open(dirname + "/particles.dat");
        for(int i=0;i<nFields ; i++)
        {
            f >> dummy;
        }

        //f.read(reinterpret_cast<char*>(_data.data()),N*getDimensions()*M*sizeof(double));
        int ii=0, tt=0;
        _data.setConstant(0);

        while(!f.eof() )
                {
                    f >> ii;
                    f >> tt;
                    for (int d=0;d<getDimensions();d++)
                    {
                        f >> _data(ii,d,tt) ;
                        _mask(tt,ii)=1;
                    }
                }
    

    f.close();
    }
    else if (format == "binary")
    {
        f.open(dirname + "/particles.dat",std::ios::binary);
        f.read(reinterpret_cast<char*>(_data.data()),N*getDimensions()*M*sizeof(double));
    }

    throw missingImplementation("Chain info not saved");
    
    // set up the groups
    particleGroups={};
    
    for (const auto & group : j["groups"])
    {
       particleGroups.push_back(
           particleGroup(
           {group["iStart"].get<int>(),group["iEnd"].get<int>(),group["iEndExtended"].get<int>(),1.0
           }
           ));
    }

    

}

void pimcConfigurations::copyData(const pimcConfigurations & confFrom, const std::array<int,2> & timeRangeFrom,const std::array<int,2> & particleRangeFrom ,
        Eigen::Tensor<Real,3> & dataTo,
       int timeOffsetTo, int particleOffsetTo )
        {
            const auto & dataFrom = confFrom.dataTensor();

            for(int t=timeRangeFrom[0], tt=timeOffsetTo;t<=timeRangeFrom[1];t++ & tt++)
            {
                for(int d=0;d<getDimensions();d++)
                   for(int i=particleRangeFrom[0], ii=particleOffsetTo;i<=particleRangeFrom[1];i++ & ii++)
                    {
                        dataTo(ii,d,tt)=dataFrom(i,d,t);
                    }
        }

        }


void pimcConfigurations::copyDataToBuffer( Eigen::Tensor<Real,2> & buffer, const std::array<int,2> & timeRange, int iParticle ,int timeOffset) const
        {

            for(int t=timeRange[0], tt=timeOffset;t<=timeRange[1];t++ & tt++)
            {
                for(int d=0;d<getDimensions();d++)
                    {
                        buffer(tt,d)=_data(iParticle,d,t);
                    }
            }
        
        }

void pimcConfigurations::copyDataFromBuffer( const Eigen::Tensor<Real,2> & buffer, const std::array<int,2> & timeRange, int iParticle ,int timeOffset)
        {

            for(int t=timeRange[0], tt=timeOffset;t<=timeRange[1];t++ & tt++)
            {
                for(int d=0;d<getDimensions();d++)
                    {
                        _data(iParticle,d,t)=buffer(tt,d);
                    }
            }
        }

void pimcConfigurations::swapData( pimcConfigurations & confFrom, const std::array<int,2> & timeRangeFrom,const std::array<int,2> & particleRangeFrom ,
        pimcConfigurations & confTo,
       int timeOffsetTo, int particleOffsetTo)
        {
             auto & dataFrom = confFrom.dataTensor();
            auto & dataTo = confTo.dataTensor();
            auto & maskFrom = confFrom.getMask();
            auto & maskTo = confTo.getMask();

            for(int t=timeRangeFrom[0], tt=timeOffsetTo;t<=timeRangeFrom[1];t++ & tt++)
            {
                for(int d=0;d<getDimensions();d++)
                   for(int i=particleRangeFrom[0], ii=particleOffsetTo;i<=particleRangeFrom[1];i++ & ii++)
                    {
                        std::swap(dataTo(ii,d,tt),dataFrom(i,d,t));
                    }
            }

        }



void pimcConfigurations::swap(int particleA, int particleB)
{
    // relink chains under the swap
    auto chainA=_chains[particleA];
    auto chainB=_chains[particleB];


    if (chainA.prev != -1)
    {
        _chains[chainA.prev].next=particleB; 
    }
    if (chainA.next != -1)
    {
        _chains[chainA.next].prev=particleB; 
    }

    if (chainB.prev != -1)
    {
        _chains[chainB.prev].next=particleA; 
    }
    if (chainA.next != -1)
    {
        _chains[chainB.next].prev=particleA; 
    }


    // swap the chain info
    std::swap(_chains[particleA],_chains[particleB]);

    // swap A and B indices in the head list
    std::replace(_heads.begin(),_heads.end(),particleA,-1);
    std::replace(_heads.begin(),_heads.end(),particleB,particleA);
    std::replace(_heads.begin(),_heads.end(),-1,particleB);


    // swap A and B indices in the tail list
    std::replace(_tails.begin(),_tails.end(),particleA,-1);
    std::replace(_tails.begin(),_tails.end(),particleB,particleA);
    std::replace(_tails.begin(),_tails.end(),-1,particleB);

    // swap data

    swapData(particleA,particleB);




}



int configurationsSampler::sampleChain(configurations_t & confs,randomGenerator_t & randG)
{
    // sample a chain with probability 1/ N_particles 
    int iParticle = uniformRealNumber(randG)*confs.nParticles();
    int k=0;
    int iChain=-1;

    for(const auto & group : confs.getGroups() )
    {
        k+=group.size();
        
        if (k> iParticle)
        {
            iChain = group.iEnd + 1 - (k-iParticle);
        }
    }

    return iChain;
}

void configurationsSampler::sampleFreeParticlePosition(
    std::array<Real,getDimensions()> & x,const std::array<Real,getDimensions()> & mean,Real tau,randomGenerator_t & randG,Real mass
){
    Real var = 2 * D * tau / mass;
    for(int d=0;d<getDimensions();d++)
    {
        x[d]=mean[d] + normal(randG)*sqrt(var);       
    }
}



};

