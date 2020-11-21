#ifndef TOOLSPIMC_H
#define TOOLSPIMC_H

#include "../src/traits.h"
#include "../src/tools.h"
#include <random>
#include <array>

namespace pimc{

    using Real = double;

    enum periodicity {periodic = 1 , open = 0};

    std::array<std::array<int,2>, 2> splitPeriodicTimeSlice(const std::array<int,2> & timeSlice, int nBeads);
    
    Real freeParticleLogProbability(std::array<Real,3> & delta,Real tau,Real mass=1);

}

#endif