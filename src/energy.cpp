#include "energy.h"
#include "wavefunction/productWavefunction.h"
#include "walkers.h"

real_t kineticEnergy::operator()(walker_t & w,wavefunction_t & psi)
	{
		real_t e=0;
		real_t ef=0;
		real_t dummy=0;

		
		psi.evaluateDerivatives(w.getStates(),w.getGradients(),w.getLogWave()     ,e,w.getTableDistances());
		
		for (const auto & grad : w.getGradients())
		{
			Eigen::Tensor<real_t,0> tmp = (grad * grad ).sum();
			ef+=tmp();	
		}
		
		return -0.5*(ef + e);
	};

real_t energy::operator()(walker_t & w,wavefunction_t & psi)
	{
		auto v=(*_pot)(w.getStates(),w.getTableDistances());

		return  kinE(w,psi) + v;
	}; 

real_t forceEnergy::operator()(walker_t & w,wavefunction_t & psi)
	{
		real_t e=0;
		real_t ef=0;
		real_t dummy=0;
		psi.evaluateDerivatives(w.getStates(),w.getGradients(),dummy,e,w.getTableDistances());
		
		for (const auto & grad : w.getGradients())
		{
			Eigen::Tensor<real_t,0> tmp = (grad * grad ).sum();
			ef+=tmp();	
		}

		auto v=(*_pot)(w.getStates(),w.getTableDistances());
		return 0.5*(ef ) + v;
	};

real_t energyFromWalker::operator()(energyFromWalker::walker_t & w,energyFromWalker::wavefunction_t & psi) 
{
  return w.getEnergy();
}
