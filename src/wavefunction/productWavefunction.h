#ifndef PRODUCT_WAVEFUNCTION_H
#define PRODUCT_WAVEFUNCTION_H

#include "wavefunction.h"
#include <initializer_list>

class productWavefunction
{
	/*
	Represent the sum of the logharitms of single wavefunctions. Does not own
	the wavefunctions.
	*/
public:
	using particles_t = wavefunction::states_t;
	using grads_t = wavefunction::grads_t;

	productWavefunction(){};

	productWavefunction(std::initializer_list<wavefunction *> waves) : _logWaves(waves) {}

	real_t operator()(const particles_t & states); // evaluates \sum_i psi_i

	void evaluateDerivatives(const particles_t & states,grads_t & grads,real_t & waveValue ,real_t & lap); // computes the gradient and computes lap: log(\prod_i psi_i) and laplacian \nabla \psi = \prod_i exp(psi_i), lap_force : sum(grad**2) 


	void evaluateDerivatives(const particles_t & states,grads_t & grads,real_t & waveValue ,real_t & lap, const tableDistances & tab); 

	void add(wavefunction * wave) {_logWaves.push_back(wave);}


	size_t size () const { return _logWaves.size();}

	const auto & waves() const {return _logWaves;}

	const wavefunction & operator()(size_t i) const {return *(_logWaves[i]);}
	wavefunction & operator()(size_t i)  {return *(_logWaves[i]);}


	const auto & getGeometry() const {return (_logWaves[0])->getGeometry();}


private:
	std::vector<wavefunction*> _logWaves;
};

#endif