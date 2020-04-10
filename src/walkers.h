#ifndef WALKERS_H
#define WALKERS_H

#include "traits.h"
#include <memory>
#include "tableDistances.h"
#include "qmcExceptions.h"
/*
A walker contains all the informiation
on the current configurations. 
Also maintains cached data for use, as particle distances , wavefunction values and
so on. Walkers own data memory
*/


struct walker
{
public:
	using grads_t = states_t;
	walker(){};
	const  auto & getStates() const {return _states;}
	const auto & getTableDistances() const {return _tab;}
	const auto & getLogWave() const {return _waveValue;}
	const auto & getGradients() const {return _gradients;}
	const auto & getLaplacianLog() const {return _lapLog;}


	auto & getStates()  {return _states;}
	auto & getTableDistances()  {return _tab;}
	auto & getLogWave() {return _waveValue;}
	auto & getGradients()  {return _gradients;}
	auto & getLaplacianLog() {return _lapLog;}
  
  virtual const real_t & getEnergy() const {throw missingImplementation("Energy not accessible from the walker"); return _waveValue;};
  virtual real_t & getEnergy()  {throw missingImplementation("Energy not accessible from the walker");return _waveValue;};
  
  
  
private:
    states_t _states; // a vector of particle data
    tableDistances _tab;
  real_t _waveValue; // value of the wavefunction
  grads_t _gradients; // contains the gradient of the wavefunction. Just a temporary
  real_t _lapLog; // contains the laplacian of the logarithm of the wavefunction
};

struct dmcWalker : public walker
{
	dmcWalker(){};
  
        virtual real_t & getEnergy() override {return _e;}
        virtual const real_t & getEnergy() const override {return _e;}
  
private:
	real_t _e=1E+20; // stores the energy of the current configuration
};

#endif
