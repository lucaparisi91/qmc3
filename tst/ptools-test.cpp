#include "ptools.h"
#include "gtest/gtest.h"
#include <random>
#include <ctime>
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
