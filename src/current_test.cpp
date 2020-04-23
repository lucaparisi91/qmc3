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
#include "dmcDriver.h"
#include "branching.h"
#include <nlohmann/json.hpp>

#include "wavefunction/jastrowWavefunctionOneBody.h"
#include "wavefunction/jastrowWavefunctionTwoBody.h"

#include "wavefunction/slaterDeterminant.h"

#include "orbitals.h"


auto getEnergies(const   pTools::walkerDistribution::walkers_t & walkers)
{
  /* Gather all energies on process 0 for this walker distribution */
  std::vector<real_t> energies;
  for (int i=0;i<walkers.size();i++)
    {
      energies.push_back(walkers[i].getEnergy() );
    }

  std::vector<int> populations;
  populations.resize(pTools::nProcesses() , 0 );
  int n = energies.size();
  MPI_Gather( & n,1 , MPI_INT,
	      populations.data(), 1, MPI_INT,
	      0, MPI_COMM_WORLD);

  
  std::vector<real_t> totalEnergies;
  int totN = 0;
  for (int i=0;i<populations.size();i++)
    {
      totN+=populations[i];
    }
  std::vector<int> offsets;
  offsets.resize(populations.size(),0);
  
  std::partial_sum(populations.begin(),populations.end() - 1,offsets.begin() + 1 );
  totalEnergies.resize(totN);
  
  MPI_Gatherv(
	      energies.data(),
	      energies.size(),
	      MPI_DOUBLE,
	      totalEnergies.data(),
	      populations.data(),
	      offsets.data(),
	      MPI_DOUBLE,
	      0,
	      MPI_COMM_WORLD
	      );
  return totalEnergies;
  
}



int main(int argc, char** argv)
{

  pTools::init(argc,argv);
  
  std::vector<int> Ns{33};
  
  real_t lBox=10000.;

  
  geometryPBC geo(lBox,lBox,lBox);
  
  real_t alpha=1.;
  
  auto J=gaussianJastrow(alpha);
  jastrowOneBodyWavefunction<gaussianJastrow> wave(J,geo,0);
	
  productWavefunction psi;
  psi.add(&wave);
  
  harmonicPotential v(geo,1.,0);
  sumPotentials pot({&v});

  energy eO(&pot);

  pTools::walkerDistribution wd;

  pTools::walkerDistribution::walkers_t walkers;

  int nP = pTools::nProcesses();
  int seed = 34;
  
  std::ranlux24 randGen(seed + pTools::rank() );
  srand(seed + pTools::rank());
  std::uniform_int_distribution<int> disWalker(7,13);

  walkers.resize( disWalker(randGen) );

  
  
  
  
  for(int i=0;i<walkers.size();i++)
    {
      state_t state(Ns[0],getDimensions() );
      state.setRandom();
      
      initializer::initialize(walkers[i],{state},psi,eO);
    }

  auto oldPopulation = wd.gatherPopulations(walkers.size());

  auto oldEnergies = getEnergies(walkers);
  
  wd.isendReceive(walkers);
  wd.wait(walkers);
  auto newPopulation = wd.gatherPopulations(walkers.size() );
  
  auto newEnergies = getEnergies(walkers);
   if (pTools::rank()==0)
     {
       int oldPopulationSize=0;
       for (auto & pop : oldPopulation)
	 {
	   oldPopulationSize+=pop;
	 }
       int newPopulationSize=0;
       for (auto & pop : newPopulation)
	 {
	   newPopulationSize+=pop;
	 }
       
       
       for (int i=0;i<newPopulation.size();i++)
	 {
	   std::cout << i << " " << oldPopulation[i] << " => " << newPopulation[i] <<  std::endl;
	 }

       std::cout << oldPopulationSize << " " << newPopulationSize << std::endl;


       for (int i=0;i<oldEnergies.size();i++)
	 {
	   std::cout << oldEnergies[i] << std::endl;
	 }
       for (int i=0;i<oldEnergies.size();i++)
	 {
	   bool found = std::find(newEnergies.begin(),oldEnergies.end(), oldEnergies[i]) != newEnergies.end() ;
	   std::cout << found << std::endl;
	 }
       
       
     }
   
  pTools::finalize();

}
