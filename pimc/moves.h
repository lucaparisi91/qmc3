#ifndef MOVESPIMC_H
#define MOVESPIMC_H


#include "../src/traits.h"
#include "pimcConfigurations.h"
#include "tools.h"
#include "metropolis.h"







namespace pimc
{

class timeSliceGenerator
{
    public:
    timeSliceGenerator(){}

    std::array<std::array<int,2>, 2> operator()(randomGenerator_t & randG, int nBeads, int maxBeadLength);
    
    private:
    std::uniform_real_distribution<float> uniformRealNumber;

};

   class levyReconstructor
    {
        public : 
        levyReconstructor(int seed, Real timeStep_) : randG(seed) , gauss(0,1), timeStep(timeStep_),_seed(seed){}

        void apply (configurations_t & configurationsNew, configurations_t & configurationsOld, int iChain , std::array<int,2> timeRange); // performs levy reconstruction on a single strand


        auto getSeed() {return _seed;}

        private:
        int _seed;
        randomGenerator_t randG;
        std::array<Real,getDimensions()> mean;
        std::normal_distribution<Real> gauss;
        Real timeStep;

    };


class firstOrderAction;


class levyMove
{

    public:

    levyMove( levyReconstructor & lev_, int maxBeadLength);

    bool attemptMove(configurations_t & confs , firstOrderAction & S);
    
    private:

    void copyToBuffer(configurations_t & confs, std::array<int,2> timeSlice,int iChain);

    void copyToBuffer(configurations_t & confs, std::array< std::array<int,2>  , 2> timeSlice,int iChain)
    {
        copyToBuffer(confs,timeSlice[0],iChain) ;
        copyToBuffer(confs,timeSlice[1],iChain); 
    };
    void copyFromBuffer(configurations_t & confs, std::array< std::array<int,2>  , 2> timeSlice,int iChain)
    {
        copyFromBuffer(confs,timeSlice[0],iChain) ;
        copyFromBuffer(confs,timeSlice[1],iChain); 
    };

    void copyFromBuffer(configurations_t & confs, std::array<int,2> timeSlice,int iChain);



    int maxBeadLength;

    levyReconstructor _levy;
    randomGenerator_t randG;
    std::uniform_real_distribution<float> uniformRealNumber;

    metropolis sampler;
    Eigen::Tensor<Real,2> tmp;
    timeSliceGenerator tGen;


};


}



#endif