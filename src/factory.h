#include "wavefunction/wavefunctionFactory.h"
#include "wavefunction/jastrowWavefunctionOneBody.h"
#include "potentialFactory.h"
#include "wavefunction/jastrowWavefunctionTwoBody.h"


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
    registerWavefunction<jastrowTwoBodyWavefunctionIndistinguishable<jastrow_t>  >() ;
  }


private:
  
  wavefunctionFactory waveFacInstance;
  potentialFactory potFacInstance;
  
};

factory&  getFactory();

std::string createId(const json_t & j);