#include "pimcConfigurations.h"
#include "qmcExceptions.h"
#include <filesystem>
#include <sstream>

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

pimcConfigurations::pimcConfigurations(
    size_t timeSlices, int dimensions, const std::vector<particleGroup> &  particleGroups_) :
    particleGroups(particleGroups_),
    M(timeSlices),N(getTotalSize(particleGroups_)),
    _data(N,dimensions,timeSlices),
     _mask( timeSlices,N),
    _nParticles(getNParticles(particleGroups_)),
    emptyGroup{0,-1,-1,0.}
{

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

int pimcConfigurations::open(   int time, int iChain ,bool copyBack)
{
    if ( _worms.size() > 0 )
    {
        throw missingImplementation("More then on worm is not supported.");

    }
    auto & group = getGroupByChain(iChain);

    if(copyBack)
    {
        for(int t=0;t<=time;t++)
        {
            for(int d=0;d<getDimensions();d++)
            {
                _data(group.iEnd+1,d,t)=_data(iChain,d,t);
            }
        }
    }
    
    group.iEnd+=1;
    int iChain2 = group.iEnd; 

    deleteBeads({0,time-1},iChain);
    createBeads({time+1,nBeads()-1},iChain);


    createBeads({0,time},iChain2);
    deleteBeads({time+1,nBeads()-1},iChain2);

    _worms.push_back({iChain2,iChain,time,time});
    return _worms.size()-1;
}

void pimcConfigurations::close( int iWorm, bool copyBack )
{
    const auto & worm = _worms[iWorm];

    auto & group = getGroupByChain(worm.iChainHead);

    // copy the data from the head to the tail
    if (copyBack)
    {
        for(int t=0;t<=worm.iTail;t++)
        {
            for(int d=0;d<getDimensions();d++)
            {
                _data(worm.iChainTail,d,t)=_data(worm.iChainHead,d,t);
            }
        }
    }
    createBeads({0,worm.iTail},worm.iChainTail  );

    
    // remove the worm from the proper group
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

void pimcConfigurations::save(const std::string & dirname)
{

    if ( ! fs::exists(dirname) ) 
    { 
        fs::create_directory(dirname); // create src folder
    }

    std::ofstream f;
    f.open(dirname + "/particles.dat",std::ios::binary);

    
    for (int t=0;t<nBeads() ;t++ )
        for (int d=0;d<getDimensions();d++)
            for (int i=0;i<nChains();i++)
            {
                f << _data(i,d,t) ;
            }
    f.close();
    nlohmann::json j;

    j["nChains"]= nChains();
    j["dimensions"]=getDimensions();
    j["timeSlices"]=nBeads();

    std::vector<nlohmann::json> wormObjects;

    for(int i=0;i<_worms.size();i++)
    {
        nlohmann::json jWorm;
        jWorm["iChain"]=_worms[i].iChainHead;
        jWorm["iHead"]=_worms[i].iHead;
        jWorm["iTail"]=_worms[i].iTail;
        wormObjects.push_back(jWorm);
    }

    j["worms"]=wormObjects;

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
    
    f.open(dirname + "/particles.dat",std::ios::binary);

    Real tmp=0;
    for (int t=0;t<nBeads() ;t++ )
        for (int d=0;d<getDimensions();d++)
            for (int i=0;i<nChains();i++)
            {
                f >> tmp ;
                _data(i,d,t) = tmp; ;
            }
    f.close();

    _mask=mask(nBeads(),nChains()) ;

    
    /*
    for (const auto & worm : j["worms"])
    {
       open( {worm["iHead"].get<int>() , worm["iTail"].get<int>() } , worm["iChain"]  );
    }
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

    */

}


};

