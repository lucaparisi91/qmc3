#include "parameters.h"
#include "wavefunction/jastrows/jastrow.h"
#include <iostream>
#include <unsupported/Eigen/CXX11/Tensor>
#include "geometry.h"
#include "wavefunction/productWavefunction.h"

using state_t = Eigen::Tensor<real_t, 2>;
using states_t = std::vector<state_t>;

real_t kinetic_energy(const real_t & lap, const states_t & grads  )
{
	real_t e=0;
	for (int i=0;i<grads.size();i++)
	{
		Eigen::Tensor<real_t,0> tmp = (grads[i] * grads[i]).sum();
		e+=tmp();
	}
	return e + lap;
}


int main(int argc, char** argv)
{
	
	int N=100;
 	state_t particleData(N , 3);
 	state_t gradient(N , 3);

 	//Eigen::Tensor<real_t, 2> diffs( (N * (N-1) )/2, 3);

 	particleData.setRandom();

 	geometryPBC geo( 10., 10., 10.);
 	auto diffs=norm(geo.differencesOneBody(particleData,{0,0,0}));

 	auto J=gaussianJastrow(1.);
 	jastrowOneBodyWavefunction<gaussianJastrow> wave(J,geo,0);

 	productWavefunction waveT;
 	waveT.add(wave);

 	real_t e , ef , waveValue =0;

 	states_t gradients {gradient};
 	states_t states {particleData};
 	waveT.evaluateDerivatives( states, gradients, waveValue, e);

 	std::cout << kinetic_energy(e, gradients) << std::endl;

}
