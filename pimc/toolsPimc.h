#ifndef TOOLSPIMC_H
#define TOOLSPIMC_H

#include "traits.h"
#include "tools.h"
#include <random>
#include <array>

namespace pimc{
#if DIMENSIONS==1
    #define TRUNCATE_D(a,b,c) a
#endif

#if DIMENSIONS==2
    #define TRUNCATE_D(a,b) a,b
#endif

#if DIMENSIONS==3
    #define TRUNCATE_D(a,b) a,b
#endif



    using json_t = nlohmann::json;
    
    using Real = double;

    enum periodicity {periodic = 1 , open = 0};

    std::array<std::array<int,2>, 2> splitPeriodicTimeSlice(const std::array<int,2> & timeSlice, int nBeads);
    
    Real freeParticleLogProbability(std::array<Real,3> & delta,Real tau,Real mass=1);

    Real average(const std::vector<Real> & observables);
    
}

#endif