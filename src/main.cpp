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
#include <string>
#include "wavefunction/jastrowWavefunctionOneBody.h"
#include "factory.h"


int main(int argc, char** argv)
{
  nlohmann::json j;
  std::cin >> j;
  std::vector<real_t> lBox;
  
  lBox=j["lBox"].get<decltype(lBox)>();
  std::vector<int> Ns;
  Ns=j["N"].get<decltype(Ns)>();
  int D=lBox.size();

  if ( D != getDimensions() )
    {
      throw invalidInput("Input file implies dimensionality different from " + std::to_string(getDimensions() ) );
      
    }
  
  geometryPBC geo( lBox[0], lBox[1], lBox[2]);
  
  states_t states;
  
  for (int i=0;i<Ns.size();i++)
    {
      state_t particleData(Ns[i] , D);
      particleData.setRandom();
      states.push_back(particleData);
    }

 
  
  getFactory().registerJastrow< gaussianJastrow >();
  
  
  auto waves = getFactory().createWavefunctions( j["wavefunctions"],geo);
  
  
  productWavefunction psi(waves);
  
  /* Potentials
     First register implemented concrete potential classes and then create the local potential as sum of individual potentials
*/

  getFactory().registerPotential<harmonicPotential>();
  getFactory().registerPotential<squareWellPotential2b>();

  auto potentials = getFactory().createPotentials(j["potentials"],geo);
  
  sumPotentials pot(potentials);
  
  energy eO(&pot);
  forceEnergy efO(&pot);
  
  realScalarEstimator m("energy",&eO);
  realScalarEstimator m2("forceEnergy",&efO);
  
  std::string method = j["method"];

  real_t timeStep = j["timeStep"];
  size_t stepsPerBlock = j["stepsPerBlock"];
  size_t nBlocks = j["nBlocks"];  
  
  if (method == "vmc")
    {
      vmcDriver vmcO(&psi,timeStep);
      vmcO.getStepsPerBlock()=stepsPerBlock;
      vmcO.getEstimators().push_back(&m);
      vmcO.getEstimators().push_back(&m2);
      vmcO.run(states,nBlocks);

    }
  else if ( method == "dmc")
    {
      size_t nW=j["walkers"];

  
      dmcDriver dmcO(&psi,&pot,timeStep,nW);
      dmcO.getStepsPerBlock()=stepsPerBlock;
      std::vector<states_t> dmcStates;
  
      for(int i=0;i<nW;i++)
	{
	  dmcStates.push_back(states);
	}
  
      dmcO.run(dmcStates,nBlocks);
    }

}
