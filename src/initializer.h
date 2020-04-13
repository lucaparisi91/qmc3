#include "traits.h"

class wavefunction;
class tableDistances;
class productWavefunction;
class walker;
class dmcWalker;
class energy;
template<class T> class walkerContainer;

struct initializer
{
  static void registerDistances(tableDistances & tab,const wavefunction & wave);
  static void registerDistances(tableDistances & tab,const productWavefunction & wave);
  
  static void initialize(walker & w, const states_t & states ,  productWavefunction & psi);
  static void initialize(dmcWalker & w, const states_t & states ,  productWavefunction & psi,energy & ob);

  static void initialize(walkerContainer<dmcWalker> & ws, const std::vector<states_t> & states ,  productWavefunction & psi,energy & ob);
  
};
