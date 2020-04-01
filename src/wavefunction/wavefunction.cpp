#include "wavefunction.h"
#include "tableDistances.h"

wavefunction::wavefunction(const geometry_t & geo_ ) : geo(&geo_)
{

};

wavefunction1b::wavefunction1b(const geometry_t & geo_ , int setA_) : wavefunction(geo_),_setA(setA_)
{
	
};


void wavefunction1b::evaluateDerivatives(const wavefunction1b::state_t & state, wavefunction1b::grad_t & gradient , real_t & wavevalue, real_t & laplacian)
{
	differences=this->getGeometry().differences(state,{0.,0.,0.});
	distances=norm(differences);

	evaluateDerivatives(state,gradient,wavevalue,laplacian,differences,distances);

}


real_t wavefunction1b::operator()(const state_t & state)
{
	differences=this->getGeometry().differences(state,{0.,0.,0.});
	distances=norm(differences);
	return (*this)(state,distances);
};

void wavefunction1b::evaluateDerivatives(const states_t & states, grads_t & gradient , real_t & wavevalue, real_t & laplacian,const tableDistances & tab)
{
	const auto & differences_from_tab = tab.differences(_setA);
	const auto & distances_from_tab = tab.distances(_setA);
	evaluateDerivatives(states[_setA], gradient[_setA],  wavevalue, laplacian ,  differences_from_tab, distances_from_tab);

}


real_t wavefunction1b::operator()(const states_t & states,tableDistances & tab ) 
{
	const auto & distances_from_tab = tab.distances(_setA);
	return (*this)(states[_setA],distances_from_tab);

}
