#ifndef MOVESPIMC_H
#define MOVESPIMC_H


#include "../src/traits.h"
#include "pimcConfigurations.h"
#include "tools.h"
#include "metropolis.h"
#include "geometryPMC.h"
#include "towerSampler.h"
#include "toolsPimc.h"

namespace pimc
{

class timeSliceGenerator
{
    public:
    timeSliceGenerator(){}

    std::array<int, 2> operator()(randomGenerator_t & randG, int nBeads, int maxBeadLength);

    private:
    std::uniform_real_distribution<float> uniformRealNumber;
};

   class levyReconstructor
    {
        public : 

        levyReconstructor( int maxReconstructionLength) :  gauss(0,1),buffer(maxReconstructionLength*2,getDimensions()) {}


        void apply (configurations_t & configurations, std::array<int,2> timeRange,int iChain ,
        Real timeStep,randomGenerator_t & randG);
        
        private:

        std::array<Real,getDimensions()> mean;
        std::normal_distribution<Real> gauss;
        Eigen::Tensor< Real , 2 > buffer;
        int _maxReconstructionLength;

    };

class firstOrderAction;

class move 
{
    public:
    virtual bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG)=0;

};


class sectorTableMoves
{
    public:
    sectorTableMoves() : totalWeight(0){}
    void push_back( move * move_,Real weight,const std::string & name="Unkown move") ; 

    bool attemptMove(configurations_t & confs, firstOrderAction & S,randomGenerator_t & randG);

    Real acceptanceRatio(int i) {return _nSuccess[i]*1./_nTrials[i];}

    std::ostream & operator>> (std::ostream & os);
    private:
    
    int sample(randomGenerator_t & random); 

    std::vector<move*> _moves;

    std::vector<Real> _nTrials;
    std::vector<Real> _nSuccess;
    std::vector<std::string> _names;

    std::vector<Real>  accumulatedWeights;
    Real totalWeight;
    towerSampler sampler;
};


class tableMoves
{
    public:

    tableMoves(){}
    
    void push_back( move * move_,Real weight,sector_t sector,const std::string & name="Unkown move") ; 

    bool attemptMove(configurations_t & confs, firstOrderAction & S,randomGenerator_t & randG);

    std::ostream & operator>> (std::ostream & os);

    private:

    sectorTableMoves openTab;
    sectorTableMoves closedTab;

    Real nOpenSectorMoves=0;
    Real nClosedSectorMoves=0;

};




class levyMove : public move
{
    public:


    levyMove( int maxBeadLength);
    levyMove(const json_t & j) : levyMove(j["reconstructionMaxLength"].get<int>() )
    {}
    bool attemptMove(configurations_t & confs , firstOrderAction & S, randomGenerator_t & randG);
    
    private:


    bool isValidSlice( configurations_t & confs , const std::array<int,2> timeRange,int iChain) const;

    int maxBeadLength;

    levyReconstructor _levy;
    std::uniform_real_distribution<float> uniformRealNumber;
    
    metropolis sampler;
    configurationsSampler confsSampler;
    Eigen::Tensor<Real,2> buffer;
    timeSliceGenerator tGen;

};


class openMove : public move
{
    public:
    // splits a chain in two morms with one overlapping bead
    openMove(Real C_ , int maxReconstructedLength_=1) ;

    openMove(const json_t & j) : openMove(j["C"].get<Real>() ,j["reconstructionMaxLength"].get<int>() ) {}


    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);
    
    private:
    Real C;
    int _maxReconstructedLength;
    std::array<Real, 3> tmp;

    const Real D = 0.5;
    configurationsSampler confsSampler;
    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    levyReconstructor _levy;
    metropolis sampler;
    Eigen::Tensor<Real,2> buffer;
};



class closeMove : public move
{
    public:
    // splits a chain in two morms with one overlapping bead

    
    closeMove(Real C_ , int maxReconstructionLength) ;
    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);

    closeMove(const json_t & j) : closeMove(j["C"].get<Real>() ,j["reconstructionMaxLength"].get<int>() ) {}


    private:
    Real C;
    std::array<Real, 3> tmp;
    int _maxLength;

    levyReconstructor _levy;
    const Real D = 0.5;
    configurationsSampler confsSampler;
    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    metropolis sampler;
    Eigen::Tensor<Real,2> buffer;
};

// advance and recede in opposite directions
class moveHead : public move
{
    public:
    moveHead(int maxAdvanceLength_);

     moveHead(const json_t & j) : moveHead(j["reconstructionMaxLength"].get<int>() ) {}

    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);

    private:

    int _maxReconstructedLength;
    std::array<Real, 3> tmp;

    const Real D = 0.5;
    configurationsSampler confsSampler;
    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    levyReconstructor _levy;
    metropolis sampler;
    Eigen::Tensor<Real,2> buffer;
};


class moveTail : public move
{
    public:
    moveTail(int maxAdvanceLength_);

    moveTail(const json_t & j) : moveTail(j["reconstructionMaxLength"].get<int>() ) {}

    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);

    private:


    int _maxReconstructedLength;
    std::array<Real, 3> tmp;

    const Real D = 0.5;
    configurationsSampler confsSampler;
    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    levyReconstructor _levy;
    metropolis sampler;
    Eigen::Tensor<Real,2> buffer;
};





// advance and recede in opposite directions
class translateMove : public move
{
    public:
    translateMove(Real max_delta, int maxBeads);

    translateMove(const json_t & j) : translateMove(j["delta"].get<Real>() ,(j["nBeads"].get<int>() + 1 ) * j["particles"][0].get<int>() ) {}
    
    
    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);

    private:

    Real _max_delta;
    std::uniform_real_distribution<Real> distr;
    Eigen::Tensor<Real,2> buffer;
    std::array<Real,getDimensions()> delta;
    configurationsSampler confSampler;
    metropolis sampler;
    std::vector<int> currentPolimerList;

};

class swapMove : public move
{
    public:
    using geometry_t = geometryPBC_PIMC;

    swapMove( int maxStepLength_, int maxN);

    swapMove( const json_t & j ) : swapMove(
        j["reconstructionMaxLength"].get<int>() ,
          j["particles"][0].get<int>()       ) {}

    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG); // will attempt to perform a new move

    private:

    int maxStepLength;
    levyReconstructor _levy;
    randomGenerator_t randG;
    std::uniform_real_distribution<float> uniformRealNumber;
    Eigen::Tensor<Real,2> buffer;
    std::vector<Real> particleSelectionAccWeights;
    Real particleSelectionWeight;
    const Real D = 0.5;
    metropolis metropolisSampler;  
    timeSliceGenerator tGen;
    configurationsSampler confSampler;
    towerSampler particleSampler;

};

template<class T>
move* __createMove(const json_t & j) {return new T(j);}


class moveConstructor
{
    public:

    moveConstructor(std::vector<int> nMaxParticles, int nBeadsMax) : _nMaxParticles(nMaxParticles),_nBeadsMax(nBeadsMax) {}

    move* createMove(const json_t & jOuter)
    {
        // creates a copy of the json input and supplement with additional essential information
        json_t j(jOuter);
        j["particles"]=_nMaxParticles;
        j["nBeads"]=_nBeadsMax;


        std::string key=j["kind"].get<std::string>(); 
        creatorMap_t::const_iterator i;
        i=creatorMap.find(key);
        if (i!= creatorMap.end())
        {
	    return (i->second)( j  );
        }
         else
        {
	    throw factoryIdNotRecorded(key);
        }
    }

    template<class T>
    void registerMove(const std::string & key)
    {
        creatorMap[key]= &__createMove<T> ; 
    }

    auto createTable(const json_t & jTable)
    {
        tableMoves tab;

        for (auto & jMove : jTable )
        {
            std::vector<std::string > sectors = jMove["sectors"];
            Real weight = jMove["weight"].get<Real>();

            auto move = createMove(jMove["move"]);

            std::string kind = jMove["move"]["kind"].get<std::string>();

            for (auto sector : sectors)
            {
                sector_t currentSector;

                if (sector == "open" )
                {
                    currentSector=sector_t::offDiagonal;
                }
                else if ( sector == "closed")
                {
                    currentSector = sector_t::diagonal;
                }
                else
                {
                    throw invalidInput("Unkown sector type");
                } 

                tab.push_back(move,weight,currentSector,kind);
            }

                 

                
            }
            

        

        return tab;

    }


    private:
    std::vector<int> _nMaxParticles;
    int _nBeadsMax;


    typedef move* (*moveCreatorFunc) ( const json_t & j);
using creatorMap_t = std::map<std::string,moveCreatorFunc>;
    creatorMap_t creatorMap;
};

//moveConstructor::creatorMap =  moveConstructor::creatorMap_t{} ;



}





#endif