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
#include "wavefunction/slaterDeterminant.h"

#include "orbitals.h"


int main(int argc, char** argv)
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

  productWavefunction psi{&wave} ;
  dmcWalker w;
  
  
  emptyPotential v(geo);
  
  energy eO(&v);
  forceEnergy efO(&v);
  
  realScalarEstimator m("energy",&eO);
  realScalarEstimator m2("forceEnergy",&efO);
  
  initializer::initialize(w,states,psi,eO);
  
  std::cout << w.getEnergy() << std::endl;

  std::cout << sineCosBasis.energy() << std::endl;
  
}
