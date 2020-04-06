#include "traits.h"

class wavefunction;
class tableDistances;
class productWavefunction;
class walker;


struct initializer
{
	static void registerDistances(tableDistances & tab,const wavefunction & wave);
	static void registerDistances(tableDistances & tab,const productWavefunction & wave);

	static void initialize(walker & w, const states_t & states ,  productWavefunction & psi);
};