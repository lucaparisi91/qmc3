
#ifndef PIMCCONFIGURATIONS
#define PIMCCONFIGURATIONS

#include <vector>
#include "../src/traits.h"
#include "../src/tools.h"
#include <unsupported/Eigen/CXX11/Tensor>
#include "qmcExceptions.h"
#include "toolsPimc.h"
#include <list>


namespace pimc
{

class mask
{
public:
    mask(){}
    mask( int nBeads, int nChains) : _mask(nChains,nBeads) {_mask.setConstant(1);}

    void setConstant(int value){_mask.setConstant(value);}

    auto & operator()(int i,int j) {return _mask(j,i);}
    const auto &  operator()(int i,int j) const {return _mask(j,i);}
private:
    Eigen::Tensor<int, 2> _mask;
};

enum sector_t{ diagonal = 0 , offDiagonal = 1 , any = 2} ;


struct particleGroup
{
    particleGroup( int iStart_,int iEnd_,int iEndExtended_,Real mass_ = 1, sector_t sector_ = sector_t::diagonal) : 
    iStart(iStart_),iEnd(iEnd_),iEndExtended(iEndExtended_),mass(mass_),sector(sector_) {}
    bool contains(int iParticle) const {return (iParticle>= iStart) and (iParticle<=iEnd);}
    int iStart; // start of the particle group
    int iEnd; // end of the active group
    int iEndExtended; // extended memory for additional particles
    Real mass;
    sector_t sector;
    auto size() const {return iEnd - iStart + 1;}  

};


class pimcConfigurations;

struct chain
{

public:
    chain() ;

    bool isOpen() const {return hasHead() or hasTail() ;}
    bool hasHead() const {return next==-1;}
    bool hasTail() const {return prev==-1;}

    bool isEmpty() const {return head -  tail <= 1; }
    void checkChainValid(); 

    int prev; // previous chain , null for a tail
    int next; // next chain, null for a head
    int head; // time Slice of the Head, will be M for closed chains and tails
    int tail; // time slice of the tail , -1 if not a a worm tail
    
};

    class pimcConfigurations
    {
        public:
        using configurationsStorage_t =  Eigen::Tensor<Real, 3> ;

        
        void saveHDF5(const std::string & filename);
        static pimcConfigurations loadHDF5(const std::string & filename);



        pimcConfigurations() : pimcConfigurations(0,getDimensions(), {} 
        ) {}


        pimcConfigurations(size_t timeSlices, int dimensions, const std::vector<particleGroup> & particleGroups_);
        bool isOpen() {return (_tails.size() > 0) or (_heads.size() > 0)  ;}


        auto &  dataTensor() {return _data;}
        const auto &  dataTensor() const {return _data;}

        auto data() {return _data.data();}
        
        auto nChains() const {return N;} 


        auto nBeads() const {return M; } 

        void fillHead(int i);

        void fillHeads();

        const auto & getChain(int i) {return _chains[i];}

        const auto & getChain(int i) const {return _chains[i];}


        const auto & getGroup(int iChain) const
        {
            for (int i=0;i<particleGroups.size();i++)
            {
                if ( particleGroups[i].contains(iChain) )
                {
                    return particleGroups[i];
                }
            }

            throw invalidInput("missing group for particle " + std::to_string(iChain));
        }


        void setHead( int iChain, int delta);
        void setTail( int iChain, int delta);


        void save(const std::string & directoryName,const std::string & format="csv") const;

        void load(const std::string & directoryName);



        int nParticles() const {return _nParticles;}

        const auto & getGroups() const {return particleGroups;}



        static void copyData(const pimcConfigurations & confFrom, const std::array<int,2> & timeRangeFrom, const std::array<int,2> & particleRangeFrom ,
          Eigen::Tensor<Real,3> & dataTo, int timeOffsetTo, int particleOffestTo );

        static void copyData(const pimcConfigurations & confFrom, const std::array<int,2> & timeRangeFrom, const std::array<int,2> & particleRangeFrom ,
         pimcConfigurations & confTo, int timeOffsetTo, int     particleOffestTo )
         {
             copyData(confFrom,timeRangeFrom,particleRangeFrom,confTo.dataTensor(),timeOffsetTo,particleOffestTo);
         }


         void copyData(const std::array<int,2> & timeRange, int iParticleFrom, int iParticleTo)
         {
             pimcConfigurations::copyData(*this,timeRange,{iParticleFrom,iParticleFrom},(*this),timeRange[0],iParticleTo);
         }

          void copyData(const std::array<int,2> & timeRange, int iParticleFrom, int timeOffset, int iParticleTo)
         {
             pimcConfigurations::copyData(*this,timeRange,{iParticleFrom,iParticleFrom},(*this),timeOffset,iParticleTo);
         }


        const auto & tails () {return _tails;}
        const auto & heads () {return _heads;}


        void swapTails(int iChainLeft,int iChainRight);

         void copyDataToBuffer(  Eigen::Tensor<Real,2> & buffer, const std::array<int,2> & timeRange, int iParticle ,int timeOffset=0) const;

         void copyDataFromBuffer(const Eigen::Tensor<Real,2> & buffer, const std::array<int,2> & timeRange, int iParticle ,int timeOffset=0);

        void swapData( pimcConfigurations & confFrom, const std::array<int,2> & timeRangeFrom,const std::array<int,2> & particleRangeFrom ,
        pimcConfigurations & confTo,
       int timeOffsetTo, int particleOffsetTo );

       void swapData(const std::array<int,2> & timeRange, int iParticleFrom, int iParticleTo)
         {
             pimcConfigurations::swapData(*this,timeRange,{iParticleFrom,iParticleFrom},(*this),timeRange[0],iParticleTo);
         }

        void swapData(int iParticleFrom, int iParticleTo)
         {
             swapData({0,nBeads()},iParticleFrom,iParticleTo);
         }

        void swap(int iParticleFrom, int iParticleTo);


        particleGroup & getGroupByChain(int iChain) {
            for (auto & group : particleGroups )
            {
                if (group.contains(iChain) )
                {
                    return  group;       
                }
            }
            throw invalidInput("Chain is not contained in any group");
        }

        
        const auto & getMask() const  {return _mask;}


        int pushChain(particleGroup & group);
        int pushChain(int iGroup);


        int next(int iChain);
        int prev(int iChain);

        int getHead(int iChain) {return _chains[iChain].head;};
        int getTail(int iChain){return _chains[iChain].tail; };

        void removeChain(int iChain); // orders of chain is not mantained

        void join( int iChainLeft, int iChainRight);

        void deleteHeadFromList(int iChain); // deletes the head from the list

        void deleteTailFromList(int iChain); // deletes the tail from the list

        std::list<int> buildPolimerList(int iChain) const; // build a list of all chains in the same permutation cycle as iChain 


        protected:


        void updateMask(const chain & chainToUpdate);

        void deleteBeads(   std::array<int,2> timeRange, int iChain ); // deactive the beads in the mask

        void createBeads(   std::array<int,2> timeRange, int iChain ); // activates the bead in the mask
        
        const auto & getChains() {return _chains; }


        int nWorms() {return (_heads.size() + _tails.size() )/2;};



        private:

        int M;
        int N;


        std::vector<particleGroup> particleGroups;
        configurationsStorage_t _data;
        std::vector<chain> _chains;

        // accelleration structures
        std::vector<int> _tails;
        std::vector<int> _heads;
        mask _mask;
        int _nParticles;
    };

    using configurations_t = pimcConfigurations; 





class configurationsSampler
{
    public:
    configurationsSampler() : uniformRealNumber(0,1),normal(0,1) {}
    
    int sampleChain(configurations_t & confs,randomGenerator_t & randG);

    void sampleFreeParticlePosition(std::array<Real,getDimensions()> & x,const std::array<Real,getDimensions()> & mean,Real tau,randomGenerator_t & randG,Real mass=1);    

    private:
    std::uniform_real_distribution<float> uniformRealNumber;
    std::normal_distribution<Real> normal;
    
    const Real D = 0.5;

};


};

#endif