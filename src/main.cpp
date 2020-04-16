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


real_t kineticEnergyGaussian(real_t alpha,distance_t dis)
{

  auto tmp=(dis* dis ).sum();

	//return -2.*alpha*alpha * tmp() +3*alpha*dis.dimensions()[0]; 

  return -2.*alpha*alpha * tmp;
}

int main(int argc, char** argv)
{
  nlohmann::json j;
  std::cin >> j;
  std::vector<real_t> lBox;
  
  lBox=j["lBox"].get<decltype(lBox)>();
  std::vector<int> Ns;
  Ns=j["N"].get<decltype(Ns)>();
  int D=lBox.size();
  
  geometryPBC geo( lBox[0], lBox[1], lBox[2]);
  
  states_t states;
  
  for (int i=0;i<Ns.size();i++)
    {
      state_t particleData(Ns[i] , D);
      particleData.setRandom();
      states.push_back(particleData);
    }
  
  tableDistances tab(geo);
  real_t alpha=1.;
  auto J=gaussianJastrow(alpha);

  jastrowOneBodyWavefunction<gaussianJastrow> wave(J,geo,0);

  productWavefunction psi{&wave};

  harmonicPotential v(geo,1.,0);
  
  energy eO(&v);
  forceEnergy efO(&v);
  
  realScalarEstimator m("energy",&eO);
  realScalarEstimator m2("forceEnergy",&efO);
	
  // vmcDriver vmcO(&psi,1e-1);
  // vmcO.getStepsPerBlock()=100000.;
  // vmcO.getEstimators().push_back(&m);
  
  // vmcO.run(states,1000);

  
  size_t nW=j["walkers"];
  real_t timeStep = j["timeStep"];
  size_t stepsPerBlock = j["stepsPerBlock"];
  size_t nBlocks = j["nBlocks"];
  dmcDriver dmcO(&psi,&v,timeStep,nW);
  dmcO.getStepsPerBlock()=stepsPerBlock;
  std::vector<states_t> dmcStates;
  
  for(int i=0;i<nW;i++)
    {
      dmcStates.push_back(states);
    }
  
  dmcO.run(dmcStates,nBlocks);	

}
