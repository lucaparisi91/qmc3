#include "traits.h"


/*
A walker contains all the informiation
on the current configurations. 
Also maintains cached data for use, as particle distances , wavefunction values and
so on. Walkers own data memory
*/



struct vmcWalker
{
	vmcWalker()
private:
    states_t _states; // a vector of particle data
	tableDistances _tab;
	real_t waveValue;
}


struct dmcWalker : vmcWalker
{
	dmcWalker()
private:
	real_t _e; // stores the energy of the current configuration
	grades_t _gradients; //  vector of the gradients of the wavefunction for each sets
};