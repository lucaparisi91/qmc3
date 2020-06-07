#include "wavefunction/wavefunctionFactory.h"
#include "wavefunction/jastrowWavefunctionOneBody.h"
#include "potentialFactory.h"
#include "wavefunction/jastrowWavefunctionTwoBody.h"
#include "wavefunction/jastrowWavefunctionThreeBody.h"

#include "wavefunction/slaterDeterminant.h"
#include "orbitals.h"
#include "estimatorsFactory.h"

/*
Defines a singleton factory  which manages the creation of wavefunctions
*/


class factory {
private:
  static factory *singleton ;
  factory(){};
  
public:
  
  static factory  * get() ;
  template<class wave_t>
  void registerWavefunction()
  {
    waveFacInstance.registerWavefunction<wave_t>();
  }
  
  template<class pot_t>
  void registerPotential()
  {
    potFacInstance.registerPotential<pot_t>();
  }

  template<class functor_t>
  void registerPotentialFunctor()
  {
    potFacInstance.registerFunctor<functor_t>();
  }
  auto createWavefunctions(json_t & j,const geometry_t & geo)
  {
    return waveFacInstance.create(j,geo);
  }

  auto createPotentials(json_t & j,const geometry_t & geo)
  {
    return potFacInstance.create(j,geo);
  }
  
  template<class jastrow_t>
  void registerJastrow()
  {
    registerWavefunction<jastrowOneBodyWavefunction<jastrow_t>  >() ;
    registerWavefunction<jastrowTwoBodyWavefunctionUndistinguishable<jastrow_t>  >() ;
    registerWavefunction<jastrowTwoBodyWavefunctionDistinguishable<jastrow_t>  >() ;
    registerWavefunction<jastrowThreeBodyWavefunctionUnDistinguishable<jastrow_t>   >() ;
    registerWavefunction<jastrowThreeBodyWavefunctionDistinguishable<jastrow_t>   >() ;

    
  }

  template<class orbital_t>
  void registerOrbital()
  {
    if (! orbital_t::isComplex() )
      registerWavefunction<slaterDeterminantWavefunction<orbitalSet<orbital_t> > >();
    else
      registerWavefunction<slaterDeterminantWavefunctionComplex<orbitalSet<orbital_t> > >();
  }

  template<class ob_t>
  void registerObservable()
  {
    estFac.registerObservable<ob_t>();
    storerFac.registerObservable<ob_t>();
    
  }

  auto createEstimators(const json_t & j)
  {
    return estFac.create(j);
  }
  
  auto createStorers(const json_t & j)
  {
    return storerFac.create(j);
  }
private:
  
  wavefunctionFactory waveFacInstance;
  potentialFactory potFacInstance;
  estimatorFactory estFac;
  storerFactory storerFac;
};

factory&  getFactory();

std::string createId(const json_t & j);
