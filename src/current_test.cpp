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


int main(int argc, char** argv)
{

  pTools::init(argc,argv);
  
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
      std::cout << w1.getStates()[0] << std::endl;
      std::cout << "-----------------" << std::endl;

      std::cout << w2.getStates()[0] << std::endl;
      
    }

  
  pTools::finalize();

}
