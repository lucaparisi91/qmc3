#include "wavefunction/productWavefunction.h"
#include "gtest/gtest.h"
#include "tableDistances.h"
#include "walkers.h"
#include "initializer.h"
#include "potential.h"
#include "energy.h"
#include "wavefunction/jastrowWavefunctionOneBody.h"
#include "wavefunction/jastrowWavefunctionTwoBody.h"

#include "wavefunction/slaterDeterminant.h"
#include "orbitals.h"
#include "estimators.h"
#include "wavefunction/jastrows/jastrowSquareWell.h"

TEST(wavefunctionTest,oneBody)
{
	int N=100;
	int D= 3;
 	state_t particleData(N , D);
 	state_t gradient(N , D);
 	real_t alpha=1.;
 	
 	particleData.setRandom();

 	geometryPBC geo( 10., 10., 10.);


 	auto J=gaussianJastrow(alpha);
 	jastrowOneBodyWavefunction<gaussianJastrow> wave(J,geo,0);
	
 	productWavefunction waveT;
 	waveT.add(&wave);
	
 	real_t e=0 , ef=0 , waveValue =0;
	
 	states_t states {particleData};
	
	walker w;
	initializer::initialize(w,states,waveT);
 	waveT.evaluateDerivatives( w);
	e=w.getLaplacianLog();
	
 	EXPECT_EQ(e,-2*alpha*N*D);

 	auto psi_value = waveT(w);
	
 	real_t sum2 = (states[0].array()*states[0].array()).sum();
 	
 	EXPECT_NEAR(psi_value,sum2*(-alpha),1e-5);

}

TEST(wavefunctionTest,tableDistances)
{
	int N=100;
	int D=3;
 	state_t particleData(N , 3);
 	state_t gradient(N , 3);


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

TEST(wavefunctionTest,2b)
{

	int N=2;
	int D=3;
 	state_t particleData(N , 3);
 	state_t gradient(N , 3);
	real_t lBox=100;
	
 	particleData.setRandom();
       
	
 	geometryPBC geo( lBox,lBox,lBox);
	
 	states_t states {particleData};

	real_t V0=2.4674011002723395 ;
	real_t R0=1.0;
	squareWellPotential2b v(geo, V0 , R0, 0., 0. ) ;

	potential_t pot({&v});
	real_t Rm=40.0;
	real_t alpha=0.5;
	
        jastrowSquareWell J(V0,R0,Rm,alpha,lBox);

	jastrowTwoBodyWavefunctionIndistinguishable<jastrowSquareWell> wave(J,geo);
	productWavefunction psi({&wave});
	dmcWalker w;

	

	energy eO(&pot);
	states[0]*=Rm/(2.*sqrt(3));
	
	initializer::initialize(w,states,psi,eO);

	double sum=0;
	const auto & dis = w.getTableDistances().distances(0,0);
	for ( int i=0;i<dis.size();i++  )
	  {
	    if (dis(i) <= R0 ) sum+=-V0;
	  }

	ASSERT_EQ(dis.size(),( N*(N-1) )/2);
	auto sum2 = v(w);
	
	EXPECT_NEAR(sum,sum2,1e-5);
	EXPECT_NEAR(eO(w,psi),0,1e-7 ) ;
	
 }



TEST(fermions,slaterWavefunctionEnergy)
{
  std::vector<int> Ns{33};

  real_t lBox=33.;
  
  geometryPBC geo(lBox,lBox,lBox);
  states_t states;
  
  for (int i=0;i<Ns.size();i++)
    {
      state_t particleData(Ns[i] , getDimensions());
      particleData.setRandom();
      particleData=particleData*lBox - lBox/2.;
      states.push_back(particleData);
    }
  
  real_t alpha=1.;
  
  orbitalSet<sinOrbital> sineCosBasis;
  
  fillFermiSea(sineCosBasis.getOrbitals(),Ns[0],lBox); // fills a fermi see with the given orbital based

  
  slaterDeterminantWavefunction<decltype(sineCosBasis)> wave(&sineCosBasis,geo,0);

  productWavefunction psi({&wave}) ;
  dmcWalker w;
  
  emptyPotential v(geo);
  sumPotentials pot({&v});
  
  energy eO(&pot);
  forceEnergy efO(&pot);
  
  realScalarEstimator m("energy",&eO);
  realScalarEstimator m2("forceEnergy",&efO);

  
  initializer::initialize(w,states,psi,eO);
  
  EXPECT_NEAR( w.getEnergy() , sineCosBasis.energy() , 1E-5 );
  
}


