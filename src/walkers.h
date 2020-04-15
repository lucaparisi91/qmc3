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



class energy;


void update(dmcWalker & w, productWavefunction & wave);
void updateForceGradientLaplacian(walker & w,productWavefunction & psi);
void updateForceGradientEnergy(dmcWalker & w,productWavefunction & psi, energy & energyOb);

template<class T>
class walkerContainer
{
public:
  using value_type=T;
  /*A wrapper around std::vector to contain walkers. 
    Calling resize does not make old allocated data invalid.
*/
  
  walkerContainer() : walkers(),_size(0) {};
  walkerContainer( std::vector<T> vec ) : walkers(vec),_size(vec.size()) {}
  
  auto & operator[](size_t i) {return walkers[i];}
  const auto & operator[](size_t i) const {return walkers[i];}
  
  void push_back( T  w)
  {
    _size=_size +1;
    
    if (_size > capacity() )
      {
	walkers.push_back(w);
      }
    else
      {
	walkers[_size-1]=w;
      }
  }
  void resize(size_t size2)
  {
    if (size2 > capacity() ) walkers.resize(size2);
    _size=size2;
  }

  void resize(size_t size2, T  w)
  {
    if (size2 > capacity() ) walkers.resize(size2,w);
    
    _size=size2;
  }
  
  size_t size() const {return _size;}
  size_t capacity() const {return walkers.size();}
  
  auto  begin()  {return walkers.begin();}
  auto end() {return walkers.begin() + _size;}

  auto  begin() const  {return walkers.begin();}
  auto end() const {return walkers.begin() + _size;}
  
  auto  cbegin() const {return walkers.cbegin();}
  auto cend() const {return walkers.cbegin() + _size;}

  
  private:
  std::vector<T> walkers;
  size_t _size;
};

#endif
