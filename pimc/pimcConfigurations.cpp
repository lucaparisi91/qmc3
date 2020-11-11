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

chain::chain(int index, int maxLength) : iChain(index) , _iTail(-1),_iHead(maxLength),_maxLength(maxLength){
    setNextChain(this);
    setPrevChain(this);
}

void chain::setHead(int i)
{
    if ( (i>_maxLength) or (i<-1) )
    {
        throw invalidInput("Head not contained on this chain");
    }
    else
    {
        _iHead=i;
    }

}

void chain::setTail(int i)
{
    if ( (i>=_maxLength) or (i<-1) )
    {
        throw invalidState("Tail not contained on this chain.");
    }
    else
    {
        _iTail=i;

    }

}

bool chain::contains( const std::array<int,2> & timeRange) const
{
    if (timeRange[0] <0 or timeRange[1]<0 )
    {
        throw invalidInput("Time indices should alwais be non negative.");
    }

    if ( getHead() > getTail()   )
    {
        if (timeRange[0] > getTail())
        {
            if ( timeRange[1] < getHead() )
            {
                return true;
            }
            else if (not hasHead() )
            {
                return getNextChain().contains({0,timeRange[1]- _maxLength});
            }

        
        }

    }
    else
    {
        if(
             (timeRange[0] < getHead() )
            and (timeRange[1] < getHead() )
        )
        {return true;} 

        if(
             (timeRange[0] > getTail() )
            and (timeRange[1] < _maxLength )
        )
            {
                return true;
            }
        else
        {
            return getNextChain().contains({0,timeRange[1]-_maxLength} );
        }
        

        


    }
    
    
    return false;
}

chain* chain::moveTail(int i)
{
if (not hasTail() )
    {
        throw invalidState("No tail to move.");
    }
    else
    {
        int newTail=getTail() + i;

        if ( hasHead() )
        {
            if ( 
                ( getHead() - newTail )*( getHead() - getTail() )<0
                 ) // ordering between chain and head changes
                 {
                    throw invalidInput("Head and tail cannot cross.");
                 }

        }
        
        
        if (newTail <= -1)
        {
            setTail(-1);

            if ( getPrevChain().hasTail() )
            {
                throw invalidState("Next chain already has a tail");
            }
            getPrevChain().setTail(_maxLength-1);
            return getPrevChain().moveTail(1 + newTail);
        }
        else
        {

        if (newTail >= _maxLength)
        {
            setTail(-1);

            if ( getNextChain().hasTail() )
            {
                throw invalidState("Next chain already has a tail");
            }

            getNextChain().setTail(0);
            return getNextChain().moveTail(newTail - _maxLength);

        }
        else
        {
            setTail(newTail);
        }
            
        }
        
    }


    return this;
}

// moves the head to the next position
void pimcConfigurations::moveHead( int iWorm, int delta )
{
    auto & currentWorm = _worms[iWorm];
    chain & oldChainHead=currentWorm.getHead();
    currentWorm.moveHead(delta);
    updateMask(oldChainHead);

    if (currentWorm.getHead().index() != oldChainHead.index() ) // head moved to an other chain
    {
        updateMask(currentWorm.getHead());
    }

};


void pimcConfigurations::moveTail( int iWorm, int delta )
{
    auto & currentWorm = _worms[iWorm];
    chain & oldChainTail=currentWorm.getTail();
    currentWorm.moveTail(delta);
    updateMask(oldChainTail);

    if (currentWorm.getTail().index() != oldChainTail.index() ) // head moved to an other chain
    {
        updateMask(currentWorm.getTail());
    }

};


chain * chain::moveHead(int i)
{
    if (not hasHead() )
    {
        throw invalidState("No head to move.");
    }
    else
    {
        int newHead=getHead() + i;

        if (newHead >= _maxLength)
        {
            setHead(_maxLength);
            if ( getNextChain().hasHead() )
            {
                throw invalidState("Next chain already has a head");
            }

            getNextChain().setHead(0);
            return getNextChain().moveHead(newHead - _maxLength);

        }
        else
        {
            if ( 
                ( newHead - getTail() )*( getHead() - getTail() )<0
                 ) // ordering between chain and head changes
                 {
                    throw invalidInput("Head and tail cannot cross.");
                 }
        }
        
    }

    return this;

}


void chain::join(chain* chainRight)
{
    auto chainLeft= this;
    chainLeft->setNextChain(chainRight);
    chainRight->setPrevChain(chainLeft);
}

worm::worm(chain * head_or_tail)
{
    if (!( head_or_tail->hasHead() xor head_or_tail->hasTail() ) )
    {
        throw invalidInput("Could not determine valid heads and tails");
    }

    if (head_or_tail->hasHead()  )
    {
        _head=head_or_tail;
        _tail=& (head_or_tail->getPrevChain());   
    }
    else 
    {
        _tail=head_or_tail;
        _head=&(head_or_tail->getNextChain());
    }

}

void worm::moveHead(int delta)
{
   _head= getHead().moveHead(delta);
}
void worm::moveTail(int delta)
{
   _tail= getTail().moveTail(delta);
}



worm::worm(chain * head,chain* tail): 
    _tail(tail),_head(head) {}


pimcConfigurations::pimcConfigurations(
    size_t timeSlices, int dimensions, const std::vector<particleGroup> &  particleGroups_) :
    particleGroups(particleGroups_),
    M(timeSlices),
    N(getTotalSize(particleGroups_)), // number of chains(includign padding chains)
    _data(N,dimensions,timeSlices),// contains a copy for buffer operations
     _mask( timeSlices,N),
    _nParticles(getNParticles(particleGroups_))// number of chains without the padding
{
    for (const auto & group : particleGroups)
    {
        for(int t=0;t<nBeads();t++)
        {
             for(int i=group.iEnd+1;i<=group.iEndExtended;i++)
             {
                 _mask(t,i)=0;
             }
            
        }
    }

    for(int i=0;i<N;i++)
    {
        _chains.emplace_back(  i, nBeads() );
    }

    
}; 

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

void pimcConfigurations::open( int time, int iChain ,bool copyBack)
{
    
    if ( _worms.size() > 0 )
    {
        throw missingImplementation("More then on worm is not supported.");
    }
    auto & group = getGroupByChain(iChain);

    group.iEnd+=1;
    int iChain2 = group.iEnd;

    // update chain information and add worm info

    auto  & chain1 = _chains[iChain];
    auto  & chain2 = _chains[iChain2];

    if (chain1.isOpen()) 
    {
        throw invalidState("The chain is already open.");
    }

    chain1.setHead(time);
    chain1.setTail(-1);

    chain2.setHead(nBeads());
    chain2.setTail(time);
    
    chain2.join(&chain1.getNextChain() );
    chain1.join(&chain2 );


    updateMask(chain1);
    updateMask(chain2);

    worm newWorm(&chain1, & chain2);
    _worms.push_back(newWorm);

    moveTail(_worms.size()-1,-1);


    if (copyBack)
    {throw missingImplementation("No copyback implemented");}

}

void pimcConfigurations::updateMask( const chain & chain)
{

    if (chain.getTail() < chain.getHead())
    {
        createBeads({chain.getTail()+1,chain.getHead()-1},chain.index() );
        deleteBeads({0,chain.getTail() },chain.index() );
        deleteBeads({chain.getHead(),nBeads() -1 },chain.index() );

    }
    else
    {
        deleteBeads({chain.getHead(),chain.getTail() },chain.index() );
        createBeads({0,chain.getHead()-1 },chain.index() );
        createBeads({chain.getTail()+1,nBeads() - 1 },chain.index() );
    }

}

void pimcConfigurations::close( int iWorm )
{
    auto & sWorm = _worms[iWorm];

    int iChainHead=sWorm.getHead().index();
    int iChainTail=sWorm.getTail().index();

    auto & group = getGroupByChain(iChainHead);

    sWorm.getHead().join(& (sWorm.getTail().getNextChain()) );
    sWorm.getHead().setHead(nBeads());
    sWorm.getHead().setTail(-1);

    updateMask(sWorm.getHead() );


    group.iEnd-=1;

    if (group.iEnd<group.iStart)
    {
        throw invalidInput("iEnd is lower than iStart. ");
    };
    
    if (_worms.size() > 1 )
    {
        throw missingImplementation("Multiple worms are not supported");
    };

    _worms.resize(0);
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


void pimcConfigurations::copyDataToBuffer( Eigen::Tensor<Real,2> & buffer, std::array<int,2> & timeRange, int iParticle ,int timeOffset) const
        {
            
            for(int t=timeRange[0], tt=timeOffset;t<=timeRange[1];t++ & tt++)
            {
                for(int d=0;d<getDimensions();d++)
                    {
                        buffer(tt,d)=_data(iParticle,d,t);
                    }
        }
        
        }


void pimcConfigurations::copyDataFromBuffer( const Eigen::Tensor<Real,2> & buffer, std::array<int,2> & timeRange, int iParticle ,int timeOffset)
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
    throw missingImplementation("Particle swap not yet supported");
}

};

