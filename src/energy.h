#include "observables.h"
#include "potential.h"
#ifndef ENERGY_H
#define ENERGY_H

class kineticEnergy : public realScalarObservable
{
public:
	virtual real_t operator()(walker_t & w,wavefunction_t & psi) override;
};

class energy : public realScalarObservable
{
public:
	energy(potential * pot_) : _pot(pot_) {}
	virtual real_t operator()(walker_t & w,wavefunction_t & psi) override;
private:
	kineticEnergy kinE;
	potential * _pot;
};

class forceEnergy : public realScalarObservable
{
	public:
	forceEnergy(potential * pot_) : _pot(pot_) {}
	virtual real_t operator()(walker_t & w,wavefunction_t & psi) override;
private:
	potential * _pot;
};

class energyFromWalker : public realScalarObservable
{
public:
  
  virtual real_t operator()(walker_t & w,wavefunction_t & psi) override;
private:
  
};

#endif
