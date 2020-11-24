#include "pimcObservables.h"
namespace pimc
{

Real thermodynamicEnergyEstimator::operator()(configurations_t & confs, firstOrderAction & S)
{
    auto & geo = S.getGeometry();

    auto & kA = S.getKineticAction();
    auto & potA = S.getPotentialAction();
    
    auto sA=kA.evaluate(confs);
    auto sV=potA.evaluate(confs);

    auto beta = confs.nBeads() * kA.getTimeStep(); 
    sA/=beta*confs.nParticles();
    sV/=beta*confs.nParticles();
    
    
    return sV - sA +  3/(2.*kA.getTimeStep());

}

};