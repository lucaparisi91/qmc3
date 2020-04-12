#include "walkers.h"
#include "energy.h"
#include "wavefunction/productWavefunction.h"

void updateForceGradientLaplacian(walker & w,productWavefunction & psi)
{
	/* Update forces ,laplacian and wavefunction value*/
  w.getTableDistances().update(w.getStates());
  psi.evaluateDerivatives(w.getStates(),w.getGradients(),w.getLogWave(),w.getLaplacianLog(),w.getTableDistances());
};

void updateForceGradientEnergy(dmcWalker & w,productWavefunction & psi, energy & energyOb)
{
  w.getTableDistances().update(w.getStates());
  w.getEnergy()=energyOb(w,psi);
};
