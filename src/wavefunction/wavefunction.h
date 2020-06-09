#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include "traits.h"
#include "geometry.h"
#include "wavefunction/jastrows/jastrow.h"
#include "qmcExceptions.h"
#include "tools.h"
#include <memory>


struct wavefunctionComponentCommands;
class tableDistances;
class walker;
class orbitalSetBase;


class constraint
{
public:
  using walker_t = walker;
  
  virtual bool operator()(const walker_t & w)=0;
  
};

class noConstraint : public constraint
{
public:
  virtual bool operator()(const walker_t & w) override  {return true;}
  
};

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
  wavefunction() {};
  virtual real_t operator()(const walker_t & state ) = 0;
  
  const auto & getGeometry() const {return *geo;}  
  
  virtual void accumulateDerivatives( walker_t & state)=0;
  
  virtual std::vector<int> sets() const = 0 ;

  virtual std::vector<orbitalSetBase*> orbitals() const {return {} ;}  
  virtual bool isComplex() const {return false;}
  
  void setGeometry(const geometry_t & geo_) {geo=&geo_;}
  
  virtual std::string print() const {return "";};

  
  virtual bool satisfyConstraints(const walker_t & state);
  
  void addConstraint(std::unique_ptr<constraint> newConstraint);
  
  
  
private:
  const geometry_t * geo;
  std::vector<std::unique_ptr<constraint> > constraints;
};




#endif
