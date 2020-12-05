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

    levyMove( levyReconstructor & lev_, int maxBeadLength);

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
    openMove(Real C_ , int maxReconstructedLength_=0) ;

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





}

#endif