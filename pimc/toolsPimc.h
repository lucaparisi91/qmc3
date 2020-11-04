#include <random>
#include "pimcConfigurations.h"

namespace pimc{

    
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