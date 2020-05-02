#include "pairCorrelation.h"
#include "tableDistances.h"
#include "walkers.h"
#include "tools.h"
#include "wavefunction/productWavefunction.h"
#include "geometry.h"
#include <iostream>

pairCorrelation::pairCorrelation(int setA_,int setB_) : setA(setA_),setB(setB_),_normalizationFactor(0)
{}

pairCorrelation::pairCorrelation(const json_t & j) : pairCorrelation(j["sets"][0],j["sets"][1])
{
  
}

void pairCorrelation::accumulate(walker_t & w,wavefunction_t & wave,pairCorrelation::accumulator_t & acc )
{
  
  auto & norms =w.getTableDistances().distances(setA,setB);
  
  setNormalizationFactor(w,wave,acc);
  for (int i=0;i<norms.rows();i++)
    {
      if (norms(i) < acc.maxx() )
  	{
  	  acc.accumulate(_normalizationFactor/(norms(i)*norms(i)),norms(i));
  	}
    }
  
  acc.weight()+=1;
}

void pairCorrelation::setNormalizationFactor(const walker_t & w , const wavefunction_t & psi ,const pairCorrelation::accumulator_t & acc) 
{
  auto lBox = psi.getGeometry().getLBox(0);
  auto dx = acc.stepSize();
  auto  NA = getN(w.getStates()[setA]);
  auto NB = getN(w.getStates()[setB]);
  _normalizationFactor=2*std::pow(lBox,3)/(dx*4*M_PI*NA*NB);
  
}
