#include "abstractFactory.h"

#include <typeinfo>

class wavefunction;

typedef wavefunction* (*wavefunctionCreatorFunc) ( const json_t & j ,geometry_t & geo);

template<class wave_t>
wavefunction * createWavefunction(const json_t & j ,geometry_t & geo)
{
  return new wave_t(j,geo);
}

class wavefunctionFactory : public abstractFactory<wavefunction,std::string, wavefunctionCreatorFunc>
{
public:
  using abstractFactory_t= abstractFactory<wavefunction,std::string, wavefunctionCreatorFunc>;

  
  template<class wave_t >
  void registerWavefunction()
  {
    registerType( wave_t::name()  , & (createWavefunction<wave_t> ) );
  }

  auto create(const json_t & j,geometry_t & geo)
  {
    std::vector<wavefunction*> waves;

    for (auto & waveJson : j )
      {
      
    
	std::string kind= waveJson["kind"];
	std::string jastrowKind = waveJson["jastrow"]["kind"];
	std::string id = kind + "/" + jastrowKind;
	waves.push_back( abstractFactory_t::create(id,waveJson,geo));
      }
    
    return waves;
  }
  

  
};

