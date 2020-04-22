#include "ptools.h"
#include "gtest/gtest.h"
#include <random>
#include <ctime>
#include "walkers.h"
#include "wavefunction/productWavefunction.h"
#include "geometry.h"
#include "traits.h"
#include "tools.h"
#include "wavefunction/jastrows/jastrow.h"
#include "wavefunction/jastrowWavefunctionOneBody.h"
#include "potential.h"
#include "energy.h"
#include "initializer.h"

TEST(pToolsTest,alias_load_balancing_algorithm_test)
{
  int Np=10;
  int Nw=94;
  int nTrials=10;
  
  int k=Nw/Np;
  
  std::vector<int> permutations;
  std::vector<int> populations(Np,k);
  
  for(int i=0;i<Nw%Np;i++)
    {
      populations[i]+=1;
    }
  
    
  std::vector<int> sources;
  std::vector<int> amounts;


  std::ranlux24  randGen(time(NULL));
  std::uniform_int_distribution<int> walkerSelector(0,Np-1);
  std::uniform_int_distribution<int> amountChange(0,int(0.5*k));
  
  
  for (int i=0;i<nTrials;i++)
    {
      int n=amountChange(randGen);
      populations[walkerSelector(randGen) ]+=n;
      do
	{
	  int i=walkerSelector(randGen);
	}
      while (populations[i]<n);
      populations[ i ]-=n;
    }
  
  std::vector<int> tmpPopulations(populations);
  std::vector<std::vector<int> > destinations;
  
  pTools::determineLoadBalanceComunicationsAliasMethod(tmpPopulations,permutations,sources,destinations,amounts);

  std::vector<int> newPopulations(populations);
  
  for ( int i =0 ; i < populations.size();i++)
    {
      // increases the current population
      newPopulations[i]+= amounts[i];
      newPopulations[  sources[i]   ]-=amounts[i];
    }
  
  auto size2 = std::accumulate(newPopulations.begin(),newPopulations.end(),0);
  
  ASSERT_EQ( Nw , size2);
  for(int i=0;i<populations.size();i++)
    {
      ASSERT_NEAR (Nw/Np, newPopulations[i] ,1);
    }
  
}

TEST(pToolsTest,walkerSend)
{
  
  if ( pTools::nProcesses() >=2 )
{
  
      std::vector<int> Ns{33};
  
  real_t lBox=10000.;
  int seed=100;

  
  geometryPBC geo(lBox,lBox,lBox);
  state_t state1(Ns[0],getDimensions() );
  state_t state2(Ns[0],getDimensions());


  srand((unsigned int) seed);

  
  state1.setRandom();
  state2.setRandom();
  
  states_t states1{state1};
  states_t states2{state2};
  
  real_t alpha=1.;
  
  auto J=gaussianJastrow(alpha);
  jastrowOneBodyWavefunction<gaussianJastrow> wave(J,geo,0);
	
  productWavefunction psi;
  psi.add(&wave);
  
  harmonicPotential v(geo,1.,0);
  sumPotentials pot({&v});

  energy eO(&pot);

  dmcWalker w1;
  dmcWalker w2;
  
  initializer::initialize(w1,states1,psi,eO);
  initializer::initialize(w2,states2,psi,eO);
  
  if (pTools::rank() == 1)
    {
      pTools::partialSend(w1,0,0);
    }
  else if (pTools::rank()==0)
    {
      pTools::partialRecv(&w2,1,0);
      
      auto & statew1 = w1.getStates()[0];
      auto & statew2= w2.getStates()[0];

      auto & gradw1 = w1.getGradients()[0];
      auto & gradw2 = w2.getGradients()[0];
      
      ASSERT_EQ(w1.getEnergy(),w2.getEnergy());
      ASSERT_EQ(w1.getLogWave(),w2.getLogWave());
      ASSERT_EQ(w1.getLaplacianLog(),w2.getLaplacianLog());

      for (int i=0;i<Ns[0];i++ )
	{
	  ASSERT_NEAR(statew1(i,0),statew2(i,0) , 1e-5);
	  ASSERT_NEAR(gradw1(i,0),gradw2(i,0) , 1e-5);
	}
      
      
    }
  
}
}


