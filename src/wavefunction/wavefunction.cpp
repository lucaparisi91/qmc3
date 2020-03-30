#include "wavefunction.h"

wavefunction::wavefunction(const geometry_t & geo_ ) : geo(&geo_)
{
	commands = new wavefunctionComponentCommands();
};

wavefunction::wavefunction(const geometry_t & geo_ , int setA) : geo(&geo_)
{
	commands = new wavefunctionSingleComponentCommands(this,setA);
};

real_t wavefunction::operator()(const wavefunction::states_t & states )
{
	return (*commands)(states);
}

void wavefunction::evaluateDerivatives(const states_t & state, grads_t & gradient , real_t & wavevalue, real_t & laplacian) 
{
	return (*commands).evaluateDerivatives(state,gradient,wavevalue,laplacian);
}

wavefunction::~wavefunction()
{
	delete commands;
}