#include "productWavefunction.h"
#include "tools.h"
#include "wavefunction.h"
#include "walkers.h"


void productWavefunction::evaluateDerivatives(productWavefunction::walker_t & w)
{
  auto & waveValue = w.getLogWave();
  auto & lap = w.getLaplacianLog();
  auto & grads = w.getGradients();
  auto & states = w.getStates();
  
  lap=0;waveValue=0;
  grads.resize(states.size());
  
  for (int i=0;i< grads.size();i++)
    {
      grads[i].resize(getN(states[i]),getDimensions());
      grads[i].setConstant(0.);
    }
  
  for (int i=0;i< size();i++)
    {
      _logWaves[i]->accumulateDerivatives(w);
    }

}



real_t productWavefunction::operator()(const productWavefunction::walker_t &states)
{
	real_t waveValue=0;
	for (const auto & wave: _logWaves) 
	  waveValue+=(*wave)(states);
	return waveValue;
}


const geometry_t & productWavefunction::getGeometry() const {return (_logWaves[0])->getGeometry();}
