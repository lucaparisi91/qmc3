#ifndef TOOLSPIMC_H
#define TOOLSPIMC_H

#include <random>
#include "pimcConfigurations.h"

namespace pimc{

    enum periodicity {periodic = 1 , open = 0};


    std::array<std::array<int,2>, 2> splitPeriodicTimeSlice(const std::array<int,2> & timeSlice, int nBeads);
    

    Real freeParticleLogProbability(std::array<Real,3> & delta,Real tau,Real mass=1);

class configurationsSampler
{
    public:
    configurationsSampler() : uniformRealNumber(0,1) {}
    
    int sampleChain(configurations_t & confs,randomGenerator_t & randG);

    void sampleFreeParticlePosition(std::array<Real,3> & x,const std::array<Real,3> & mean,Real tau,randomGenerator_t & randG,Real mass=1);
    

    private:
    std::uniform_real_distribution<float> uniformRealNumber;
    std::normal_distribution<Real> normal;

    const Real D = 0.5;

};



}

#endif