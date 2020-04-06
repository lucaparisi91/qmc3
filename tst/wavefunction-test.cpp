#include "wavefunction/productWavefunction.h"
#include "gtest/gtest.h"
#include "tableDistances.h"
#include "walkers.h"
#include "initializer.h"
#include "potential.h"
#include "energy.h"


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


 	auto J=gaussianJastrow(alpha);
 	jastrowOneBodyWavefunction<gaussianJastrow> wave(J,geo,0);

 	productWavefunction waveT;
 	waveT.add(&wave);

 	real_t e , ef , waveValue =0;

 	states_t gradients {gradient};
 	states_t states {particleData};
 	waveT.evaluateDerivatives( states, gradients, waveValue, e);
 	
 	EXPECT_EQ(e,-2*alpha*N*D);

 	auto psi_value = waveT(states);

 	Eigen::Tensor<real_t,0> sum2 = (states[0]*states[0]).sum();
 	
 	EXPECT_NEAR(psi_value,sum2()*(-alpha),1e-5);

}

TEST(wavefunctionTest,tableDistances)
{
	int N=100;
	int D=3;
 	state_t particleData(N , 3);
 	state_t gradient(N , 3);

 	//Eigen::Tensor<real_t, 2> diffs( (N * (N-1) )/2, 3);

 	particleData.setRandom();

 	geometryPBC geo( 10., 10., 10.);

 	states_t states {particleData,2*particleData};
 	tableDistances tab(geo);

 	tab.add(0);

 	tab.add(0,0);
 	tab.add(0,1);
 	tab.update(states);

 	const auto & differences = tab.differences(0,0);

 	int k=0;
 	for (int i=0;i<N;i++)
 	{
 		for(int j=0;j<i;j++)
 		{
 			for(int d=0;d<D;d++)
 			{
 				auto diff1=geo.difference(particleData(i,d) - particleData(j,d),d ) ;
 				auto diff2 = differences(k,d);

 				EXPECT_NEAR(diff1, diff2, 1e-4);	

 			}
 				k++;
 		}

 	}

}

TEST(wavefunctionTest,oneBody_from_distances)
{
	int N=100;
	int D= 3;
 	state_t particleData(N , D);
 	state_t gradient(N , D);
 	real_t alpha=1.;
 	//Eigen::Tensor<real_t, 2> diffs( (N * (N-1) )/2, 3);

 	particleData.setRandom();

 	geometryPBC geo( 10., 10., 10.);


 	auto J=gaussianJastrow(alpha);
 	jastrowOneBodyWavefunction<gaussianJastrow> wave(J,geo,0);

 	productWavefunction waveT;
 	waveT.add(&wave);

 	real_t e , ef , waveValue =0;

 	states_t gradients {gradient};
 	states_t states {particleData};

 	tableDistances tab(geo);

 	tab.add(0);

 	tab.update(states);

 	waveT.evaluateDerivatives( states, gradients, waveValue, e,tab);
 	
 	EXPECT_EQ(e,-2*alpha*N*D);
 	
}

TEST(wavefunctionTest,harmonic_oscillator_3d)
{

	int N=100;
	int D=3;
 	state_t particleData(N , 3);
 	state_t gradient(N , 3);

 	//Eigen::Tensor<real_t, 2> diffs( (N * (N-1) )/2, 3);

 	particleData.setRandom();

 	geometryPBC geo( 10., 10., 10.);

 	states_t states {particleData};
 	tableDistances tab(geo);
 	real_t alpha=0.5;
 	auto J=gaussianJastrow(alpha);

 	jastrowOneBodyWavefunction<gaussianJastrow> wave(J,geo,0);

 	productWavefunction psi{&wave};

 	walker w;
 	initializer::initialize(w,states,psi);

 	harmonicPotential v(geo,1.,0);

 	energy eO(&v);

 	auto e = eO(w,psi);
 	//auto ek = kineticEnergyGaussian(alpha, w.getTableDistances().distances(0));
 	EXPECT_NEAR(e,150.,1e-5);
 	
 }