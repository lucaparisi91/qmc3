#include "parameters.h"
#include "wavefunction/jastrows/jastrow.h"
#include <iostream>
#include <unsupported/Eigen/CXX11/Tensor>
#include "geometry.h"
#include "wavefunction/productWavefunction.h"
#include "tableDistances.h"
#include "potential.h"
#include "initializer.h"
#include "walkers.h"
#include "energy.h"
#include "estimators.h"
#include "driver.h"

using state_t = Eigen::Tensor<real_t, 2>;
using states_t = std::vector<state_t>;


real_t kineticEnergyGaussian(real_t alpha,distance_t dis)
{

	Eigen::Tensor<real_t,0> tmp=(dis * dis).sum();

	//return -2.*alpha*alpha * tmp() +3*alpha*dis.dimensions()[0]; 

	return -2.*alpha*alpha * tmp();
}


int main(int argc, char** argv)
{
	int N=1;
	int D=3;
 	state_t particleData(N , 3);
 	state_t gradient(N , 3);

 	//Eigen::Tensor<real_t, 2> diffs( (N * (N-1) )/2, 3);

 	particleData.setRandom();

 	geometryPBC geo( 100., 100., 100.);

 	states_t states {particleData};
 	tableDistances tab(geo);
 	real_t alpha=2.;
 	auto J=gaussianJastrow(alpha);

 	jastrowOneBodyWavefunction<gaussianJastrow> wave(J,geo,0);

 	productWavefunction psi{&wave};

 	walker w;
 	initializer::initialize(w,states,psi);

 	harmonicPotential v(geo,1.,0);

 	energy eO(&v);
 	forceEnergy efO(&v);


 	realScalarEstimator m("energy",&eO);
 	realScalarEstimator m2("forceEnergy",&efO);

 	

 	auto e = eO(w,psi);
 	//auto ek = kineticEnergyGaussian(alpha, w.getTableDistances().distances(0));
 	std::cout << e << std::endl;
 	//std::cout << ek << std::endl;

 	vmc::vmcDriver vmcO(&psi,1e-2);
 	vmcO.stepsPerBlock()=100000.;
 	vmcO.estimators().push_back(&m);
 	vmcO.estimators().push_back(&m2);

 	vmcO.run(states,10000);

}
