#include "wavefunction/productWavefunction.h"
#include "gtest/gtest.h"

using state_t = Eigen::Tensor<real_t, 2>;
using states_t = std::vector<state_t>;

TEST(wavefunctionTest,oneBody)
{
	int N=100;
	int D= 3;
 	state_t particleData(N , D);
 	state_t gradient(N , D);
 	real_t alpha=1.;
 	//Eigen::Tensor<real_t, 2> diffs( (N * (N-1) )/2, 3);

 	particleData.setRandom();

 	geometryPBC geo( 10., 10., 10.);
 	auto diffs=norm(geo.differencesOneBody(particleData,{0,0,0}));

 	auto J=gaussianJastrow(alpha);
 	jastrowOneBodyWavefunction<gaussianJastrow> wave(J,geo,0);

 	productWavefunction waveT;
 	waveT.add(wave);

 	real_t e , ef , waveValue =0;

 	states_t gradients {gradient};
 	states_t states {particleData};
 	waveT.evaluateDerivatives( states, gradients, waveValue, e);
 	
 	EXPECT_EQ(e,-2*alpha*N*D);

}