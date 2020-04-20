#include "wavefunction.h"

template<class orbitalSet_t>
class slaterDeterminantWavefunction : public wavefunction
{
public:
  using grad_t = ::difference_t;
  
  slaterDeterminantWavefunction(orbitalSet_t * orbital,geometry_t & geo,size_t setA);

  virtual void accumulateDerivatives(walker_t & w);

  virtual real_t operator()(const walker_t & w);

  virtual std::vector<int> sets() const override {return {setA};}
  
  std::vector<orbitalSetBase*> orbitals() const override {return {_orbitals};}

  static std::string name() {return "slater/" + orbitalSet_t::name();}

private:
  grad_t tmpGrad;
  orbitalSet_t *_orbitals;
  int setA;
};
