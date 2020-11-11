
#ifndef PIMCCONFIGURATIONS
#define PIMCCONFIGURATIONS

#include <vector>
#include "../src/traits.h"
#include "../src/tools.h"
#include <unsupported/Eigen/CXX11/Tensor>
#include "qmcExceptions.h"

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

struct particleGroup
{
    
    bool contains(int iParticle) const {return (iParticle>= iStart) and (iParticle<=iEnd);}
    int iStart; // start of the particle group
    int iEnd; // end of the active group
    int iEndExtended; // extended memory for additional particles
    Real mass;

    auto size() const {return iEnd - iStart + 1;}
    
};


struct chain
{

public:

    chain(int index, int max_length);
        
    void join(chain * chainRight);
    
    bool isOpen() const {return hasHead() or hasTail() ;}
    bool hasHead() const {return _iHead != _maxLength;}
    bool hasTail() const {return _iTail!=-1;}

    bool contains(const std::array<int,2> & timeRange) const ;

    


    chain &  getNextChain()  {return *_nextChain; }
    chain &  getPrevChain() {return *_prevChain; }

    const chain &  getNextChain() const  {return *_nextChain; }
    const chain &  getPrevChain() const {return *_prevChain; }


    void setNextChain(chain * chainNew)  {_nextChain=chainNew;}
    void setPrevChain(chain * chainNew) {_prevChain=chainNew;}


    void setHead (int i);
    void setTail(int i);


    chain * moveHead(int delta); // returns the location of the new head 
    chain * moveTail(int delta); // retruns the location of the new tail

    
    int  index() const {return iChain; }

    const int & getHead () const {return _iHead; }
    const int & getTail  () const {return _iTail; }

    private:    

    int iChain;
    chain* _prevChain; // previous chain , null for a tail
    chain* _nextChain; // next chain, null for a head
    int _iHead; // time Slice of the Head, will be M for closed chains and tails
    int _iTail; // time slice of the tail , -1 if not a a worm tail
    int _maxLength; // maximum number of beads

};

struct worm
{
    worm() : worm(nullptr,nullptr) {}
    worm(chain* head_or_tail);

    worm(chain * head, chain * tail) ;

    auto &  getTail() {return *_tail;}
    auto &  getHead() {return *_head;}

    const auto &  getTail() const {return *_tail;}
    const auto &  getHead() const {return *_head;}


    void moveHead(int delta); // moves the head by delta steps (positive or negative)
    void moveTail(int delta);

    private:

    chain * _tail;
    chain * _head;
};

    class pimcConfigurations
    {
        public:
        using configurationsStorage_t =  Eigen::Tensor<Real, 3> ;

        pimcConfigurations() : pimcConfigurations(0,getDimensions(), {} 
        ) {}

        pimcConfigurations(size_t timeSlices, int dimensions, const std::vector<particleGroup> & particleGroups_);
        bool isOpen() {return _worms.size()> 0 ;}

        auto &  dataTensor() {return _data;}
         const auto &  dataTensor() const {return _data;}


        auto data() {return _data.data();}
        
        auto nChains() const {return N;} 
        
        auto nBeads() const {return M; } 


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

        void open( int time, int iChain,bool copyBack=false); // adds a worm
        void close(int iWorm); // closes the ith worm (not the the ith chain)


        const auto & worms() const {return _worms;} // do not allows worms to be manipulated outside of class to keep the mask in sync


        void moveHead( int iWorm, int delta);
        void moveTail( int iWorm, int delta);


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


         void copyDataToBuffer( Eigen::Tensor<Real,2> & buffer, std::array<int,2> & timeRange, int iParticle ,int timeOffset=0) const;


         void copyDataFromBuffer(const  Eigen::Tensor<Real,2> & buffer, std::array<int,2> & timeRange, int iParticle ,int timeOffset=0);

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

        void setHead(chain & chain, int time);


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

        void join( int chainLeft, int chainRight)
        {
            _chains[chainLeft].join(&( _chains[chainRight]));
        }


        const auto & getMask() const  {return _mask;}

        const auto & getChainsInfo() const {return _chains;}

        protected:
        void updateMask(const chain & chainToUpdate);


        void deleteBeads(   std::array<int,2> timeRange, int iChain ); // deactive the beads in the mask

        void createBeads(   std::array<int,2> timeRange, int iChain ); // activates the bead in the mask


        private:
        


        int M;
        int N;

        std::vector<particleGroup> particleGroups;
        std::vector<worm> _worms; // order is not preserved during closes
        configurationsStorage_t _data;
        std::vector<chain> _chains;
        mask _mask;
        int _nParticles;

    };

    using configurations_t = pimcConfigurations; 

};

#endif