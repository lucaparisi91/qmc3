#include "slaterDeterminant.h"
#include "slaters.h"
#include "orbitals.h"
#include "walkers.h"


template<class orbitalSet_t>
slaterDeterminantWavefunction<orbitalSet_t>::slaterDeterminantWavefunction(orbitalSet_t * orbitals_,geometry_t & geo_,size_t setA_) : _orbitals(orbitals_),setA(setA_),wavefunction::wavefunction(geo_) {}


template<class orbitalSet_t>
real_t  slaterDeterminantWavefunction<orbitalSet_t>::operator()(const walker_t & w)
{
  return w.getTableSlaters().logDeterminant(setA);
};


template<class orbitalSet_t>
void slaterDeterminantWavefunction<orbitalSet_t>::accumulateDerivatives(walker_t & w)
{
 
  const auto & inverseSlater = w.getTableSlaters().slaterMatrixInverse(setA);
  const auto & state = w.getStates()[setA];
  const auto & orbitals= _orbitals->getOrbitals();
  auto & lap=w.getLaplacianLog();
  
  auto & grad = w.getGradients()[setA];
  
  const auto & slater = w.getTableSlaters().slaterMatrix(setA);
  real_t tmp2=0;
  
  const int N = getN(state);
  
  constexpr int D = getDimensions();
  std::array<real_t,3> dx;
  assert( N <= orbitals.size());

  tmpGrad.resize(N,D);
  tmpGrad.setConstant(0);
  for(int j=0;j<N;j++)
      {
	const auto & orbital = orbitals[j];
	
	for(size_t i=0;i<N;i++)
	  {
	    
	    orbital.evaluateDerivatives( state(i,0),state(i,1),state(i,2) , dx[0],dx[1],dx[2],tmp2);
	    for(size_t id=0;id<D;id++)
	      {
		tmpGrad(i,id)+=dx[id]*inverseSlater(j,i);
	      }
	    
	    lap+=tmp2*inverseSlater(j,i);
	    
		
	  }
	    

	    
		    	  
      }
  lap-= (tmpGrad*tmpGrad).sum();
  grad+=tmpGrad;
  w.getLogWave()+=(*this)(w);
}


template class slaterDeterminantWavefunction<orbitalSet<sinOrbital> > ;
