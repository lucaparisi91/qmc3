#include "trimerCorrelation.h"
#include "walkers.h"
#include "wavefunction/productWavefunction.h"
#include "geometry.h"
#include "tools.h"

void trimerCorrelation::accumulate(walker_t & w,wavefunction_t & wavefunction,accumulator_t & acc)
{
  assert(getN(w.getStates()[setC] )== 1);

  setNormalizationFactor(w,wavefunction,acc);
    auto & norms1 =w.getTableDistances().distances(setA,setC);
    auto & norms2 =w.getTableDistances().distances(setB,setC);

    auto & norms3 = w.getTableDistances().distances(setA,setB);
    real_t radius=0;

    auto N1 = norms1.size();
    auto N2 = norms2.size();
    auto N3 = norms2.size();


     
  for (int i=0;i<N1;i++ )
    for (int j=0;j<N2;j++ )
      {
	auto d3 =norms3(i*N2 + j);
	radius=std::sqrt( norms1(i)*norms1(i) + norms2(j)*norms2(j) + d3*d3 );
	if ( radius < acc.maxx() )
	  {
#if DIMENSIONS == 1
	      acc.accumulate(_normalizationFactor,radius);
#endif
#if DIMENSIONS == 3
	      
	    acc.accumulate(_normalizationFactor/(radius*radius),radius);
#endif
	  }
      }
  
  acc.weight()+=1;
}


trimerCorrelation::trimerCorrelation(int setA_,int setB_,int setC_) : setA(setA_),setB(setB_),setC(setC_),_normalizationFactor(0)
{
    assert(setA!=setB);
    assert(setA!=setC);
}

trimerCorrelation::trimerCorrelation(const json_t & j) : trimerCorrelation(j["sets"][0],j["sets"][1],j["sets"][2])
{
  
}



void trimerCorrelation::setNormalizationFactor(const walker_t & w , const wavefunction_t & psi ,const trimerCorrelation::accumulator_t & acc) 
{
  auto lBox = psi.getGeometry().getLBox(0);
  auto dx = acc.stepSize();
  auto  NA = getN(w.getStates()[setA]);
  auto NB = getN(w.getStates()[setB]);
  auto NC = getN(w.getStates()[setB]);
  
#if DIMENSIONS == 3
  _normalizationFactor=1/(dx*4*M_PI*NA*NB*NC);
#endif

#if DIMENSIONS == 1
  
  _normalizationFactor=std::pow(lBox,1)/(dx*NA*NB*NC*2);   
#endif
  
  
}
