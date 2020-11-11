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
        levyReconstructor( int maxReconstructionLength) :  gauss(0,1),buffer(maxReconstructionLength,getDimensions()) {}
        

        void apply (configurations_t & configurationsNew, configurations_t & configurationsOld , std::array<int,2> timeRange, int iChain,Real timeStep,randomGenerator_t & randG); // performs levy reconstruction on a single strand


        void apply (configurations_t & configurationsNew, std::array<int,2> timeRange,int iChain ,
        Real timeStep,randomGenerator_t & randG)
        {
            apply(configurationsNew,configurationsNew,timeRange,iChain,timeStep,randG);
        }
        
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


class tableMoves
{
    public:
    tableMoves() : totalWeight(0){}
    void push_back( move * move_,Real weight) ; 

    move & sample(randomGenerator_t & random); // selects a move

    private:
    std::vector<move*> _moves;
    std::vector<Real>  accumulatedWeights;
    Real totalWeight;
    towerSampler sampler;
};

class levyMove : public move
{
    public:

    levyMove( levyReconstructor & lev_, int maxBeadLength);

    bool attemptMove(configurations_t & confs , firstOrderAction & S, randomGenerator_t & randG);
    
    private:


    bool isValidSlice( configurations_t & confs , const std::array<int,2> timeRange,int iChain) const;
    void copyToBuffer(configurations_t & confs, std::array<int,2> timeSlice,int iChain,int offset=0);

    void copyToBuffer(configurations_t & confs, std::array< std::array<int,2>  , 2> timeSlice,int iChain)
    {
        copyToBuffer(confs,timeSlice[0],iChain) ;
        copyToBuffer(confs,timeSlice[1],iChain,timeSlice[0][1]-timeSlice[0][0] + 1); 
    };

    void copyFromBuffer(configurations_t & confs, std::array< std::array<int,2>  , 2> timeSlice,int iChain)
    {
        copyFromBuffer(confs,timeSlice[0],iChain) ;
        copyFromBuffer(confs,timeSlice[1],iChain, timeSlice[0][1]-timeSlice[0][0] + 1);

    };


    void copyFromBuffer(configurations_t & confs, std::array<int,2> timeSlice,int iChain,int offset=0);


    int maxBeadLength;

    levyReconstructor _levy;
    std::uniform_real_distribution<float> uniformRealNumber;
    
    metropolis sampler;
    configurationsSampler confsSampler;
    Eigen::Tensor<Real,2> tmp;
    timeSliceGenerator tGen;

};


class openMove
{
    public:
    // splits a chain in two morms with one overlapping bead
    openMove(Real C_ , Real timeStep_) : C(C_),
    timeStep(timeStep_) {}

    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);


    private:
    Real C;
    int maxBeadLength;
    std::array<Real, 3> tmp;
    Real timeStep;

    const Real D = 0.5;
    configurationsSampler confsSampler;
    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;

    metropolis sampler;
};

class closeMove
{
    public:
    // splits a chain in two morms with one overlapping bead
    closeMove(Real C_ , Real timeStep_) : C(C_),
    timeStep(timeStep_) {}

    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);
    
    private:
    Real C;
    int maxBeadLength;
    std::array<Real, 3> tmp;
    Real timeStep;

    const Real D = 0.5;
    configurationsSampler confsSampler;
    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    metropolis sampler;
};


// advance and recede in opposite directions
class advanceRecedeMove
{
    public:
    advanceRecedeMove(int maxAdvanceLength_);

    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);
    private:

    int maxAdvanceLength;
    timeSliceGenerator tGen;
    levyReconstructor _levy;
    std::uniform_real_distribution<Real> distr;
    Eigen::Tensor<Real,2> buffer;

    configurationsSampler confSampler;
    metropolis sampler;
    std::array<Real,3> tmpPosition;
    std::array<Real,3> tmpMean;

};


class swapMove
{
    public:
    using geometry_t = geometryPBC_PIMC;

    swapMove(Real rmax_, int maxStepLength_, int seed);
    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG); // will attempt to perform a new move


    private:
    Real rmax;
    Real sigma;

    int maxStepLength;
    levyReconstructor _levy;
    randomGenerator_t randG;
    std::uniform_real_distribution<float> uniformRealNumber;
    geometry_t geo;
    Eigen::Tensor<Real,2> tmpChainWorm;
    Eigen::Tensor<Real,2> tmpChainParticle;
    std::vector<Real> particleSelectionAccWeights;
    Real particleSelectionWeight;
    const Real D = 0.5;
    towerSampler particleSampler;
    metropolis metropolisSampler;  
    timeSliceGenerator tGen;
    configurationsSampler confSampler;
    towerSampler chainSampler;

};



}

#endif