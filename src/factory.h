#include "wavefunction/wavefunctionFactory.h"
#include "wavefunction/jastrowWavefunctionOneBody.h"

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

  auto createWavefunctions(json_t & j,geometry_t & geo)
  {
    return waveFacInstance.create(j,geo);
  }

  
  template<class jastrow_t>
  void registerJastrow()
  {
    registerWavefunction<jastrowOneBodyWavefunction<jastrow_t>  >() ;
  }


private:
  
  wavefunctionFactory waveFacInstance;
  
};

factory&  getFactory();
