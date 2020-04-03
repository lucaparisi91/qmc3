#include "traits.h"
#include <memory>
#include "tableDistances.h"

/*
A walker contains all the informiation
on the current configurations. 
Also maintains cached data for use, as particle distances , wavefunction values and
so on. Walkers own data memory
*/


struct walker
{
	using grads_t = states_t;
	walker(){};
	const  auto & getStates() const {return _states;}
	const auto & getTableDistances() const {return _tab;}
	const auto & getLogWave() const {return _waveValue;}
	const auto & getGradients() const {return _gradients;}

	auto & getStates()  {return _states;}
	auto & getTableDistances()  {return _tab;}
	auto & getLogWave() {return _waveValue;}
	auto & getGradients()  {return _gradients;}

private:
    states_t _states; // a vector of particle data
	tableDistances _tab;
	real_t _waveValue; // value of the wavefunction
	grads_t _gradients; // contains the gradient of the wavefunction
};


struct dmcWalker : walker
{
	dmcWalker();
	const auto & energy() {return _e;}
private:
	real_t _e=1E+20; // stores the energy of the current configuration
};