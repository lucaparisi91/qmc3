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

  auto initialConfigurations = readStatesFromDirectory("configurations");
  
  pTools::finalize();

}
