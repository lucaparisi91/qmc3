#include "potential.h"
#include "geometry.h"
#include "tableDistances.h"
#include <iostream>
potential::potential(const geometry_t & geo_ ) : geo(&geo_)
{

};

potential1b::potential1b(const geometry_t & geo_ , int setA_) : potential::potential(geo_),_setA(setA_)
{
	
};

real_t potential1b::operator()(const state_t & state)
{
	differences=this->getGeometry().differences(state,{0.,0.,0.});
	distances=norm(differences);
	return (*this)(state,distances);
};

real_t potential1b::operator()(const states_t & states,const tableDistances & tab ) 
{
	const auto & distances_from_tab = tab.distances(_setA);
	return (*this)(states[_setA],distances_from_tab);
}


harmonicPotential::harmonicPotential(const geometry_t & geo,real_t freq , int setA ) : potential1b(geo,setA) ,omega(freq) {}

real_t harmonicPotential::operator()(const state_t & state,const distance_t & dis)
{
		
	
		int N=state.dimensions()[0];
		real_t sum=0;

		for(int i=0;i<N;i++)
		{
			sum+=0.5 * omega * omega * dis(i) * dis(i);
		}
		return sum;
};