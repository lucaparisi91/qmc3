#include "productWavefunction.h"


void productWavefunction::evaluateDerivatives(const productWavefunction::particles_t & states, productWavefunction::grads_t & grads, real_t & waveValue,real_t & lap)
{
	waveValue=0;lap=0;

	real_t lapPartial,waveValuePartial;
	for (int i=0;i< size();i++)
	 {
	 	grads[i].setConstant(0.);
	 }
	
	for (int i=0;i< size();i++)
		{
			_logWaves[i]->evaluateDerivatives(states,grads,waveValuePartial,lapPartial);
			waveValue+=waveValuePartial;
			lap+=lapPartial;
		}

}

real_t productWavefunction::operator()(const productWavefunction::particles_t &states)
{
	real_t waveValue=0;
	for (const auto & wave: _logWaves) 
		 waveValue+=(*wave)(states);
	return waveValue;
}