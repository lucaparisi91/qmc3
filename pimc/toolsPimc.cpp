#include "toolsPimc.h"

namespace pimc
{

int configurationsSampler::sampleChain(configurations_t & confs,randomGenerator_t & randG)
{
    // sample a chain with probability 1/ N_particles 
    int iParticle = uniformRealNumber(randG)*confs.nParticles();
    int k=0;
    int iChain=-1;

    for(const auto & group : confs.getGroups() )
    {
        k+=group.size();
        
        if (k> iParticle)
        {
            iChain = group.iEnd + 1 - (k-iParticle);
        }
    }
    
    return iChain;
}

void configurationsSampler::sampleFreeParticlePosition(
    std::array<Real,getDimensions()> & x,const std::array<Real,getDimensions()> & mean,Real tau,randomGenerator_t & randG,Real mass
){
    Real var = 2 * D * tau / mass;
    for(int d=0;d<getDimensions();d++)
    {
        x[d]=mean[d] + normal(randG)*sqrt(var);       
    }
}

Real freeParticleLogProbability(std::array<Real,3> & delta,Real tau,Real mass)
    {
        const Real D = 0.5;
        Real var = 2 * D * tau / mass;

        Real p=0;
        
        for(int d=0;d<getDimensions();d++)
        {
            p+= delta[d]*delta[d];
        }
        p*=-0.5 /var;
        p+= -0.5*log(2*M_PI*var);
        
        return p;
    }




}