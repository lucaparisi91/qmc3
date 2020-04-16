#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include "traits.h"
#include "geometry.h"
#include "wavefunction/jastrows/jastrow.h"
#include "qmcExceptions.h"
#include "tools.h"

struct wavefunctionComponentCommands;
class tableDistances;
class walker;
class orbitalSetBase;

class wavefunction
{
	/*
	Represents the logarithm of the wavefunction. Supports wavefunctions acting on at most two different sets
	*/
public:
  using geometry_t = geometry;
  using state_t = ::state_t;
  using grad_t = state_t;
  using states_t= ::states_t;
  using grads_t = ::states_t;
  using walker_t = walker;
  
  wavefunction(const geometry_t & geo_ );
  
  virtual real_t operator()(const walker_t & state ) = 0;
  
  const auto & getGeometry() const {return *geo;}  
  
  virtual void accumulateDerivatives( walker_t & state)=0;
  
  virtual std::vector<int> sets() const = 0 ;

  virtual std::vector<orbitalSetBase*> orbitals() const {return {} ;}
  
private:
  const geometry_t * geo;

};




#endif
