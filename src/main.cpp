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
#include "dmcDriver.h"
#include "vmcDriver.h"
#include "branching.h"
#include <string>
#include "wavefunction/jastrowWavefunctionOneBody.h"
#include "factory.h"
#include "wavefunction/jastrows/jastrowSquareWell.h"
#include "ptools.h"
#include "pairCorrelation.h"
#include "centerOfMassSquared.h"
#include "hyperRadius.h"


bool check_n_particles(const states_t & states,const std::vector<int> & Ns)
{
  bool pass=true;
  for (int i=0;i<Ns.size();i++)
    {
      pass=pass and ( getN(states[i]) == Ns[i] );
    }
  return pass;
}


template<class T>
auto   find( json_t & j,std::string key,const T & value)
{
  for (auto it = j.begin() ; it<j.end() ; it++ )
    {
      if ( (*it)[key] == value )
	{
	  return it;
	}
    }
  
  return j.end();
}


int main(int argc, char** argv)
{

  /*
    Reads input from standrard input and broadcast to all other mpi processes
*/
  pTools::init(argc,argv);

  std::string jSonString="";  
  
  
  if ( pTools::isMaster())
    {
      std::string line;
      while (std::getline(std::cin, line))
	{
	  jSonString+=line;
	}
  
    }
  
  pTools::broadcast(&jSonString,0);
  
  std::istringstream ss(jSonString);
  
  nlohmann::json j;

  ss >> j;

  
  

  if (pTools::isMaster())
     {
       std::cout << "MPI processes: " <<  pTools::nProcesses() << std::endl;
     }
  
  
  std::vector<real_t> lBox;
  std::vector<real_t> lBoxInitialCondition;
  lBox=j["lBox"].get<decltype(lBox)>();
  std::vector<int> Ns;
  Ns=j["N"].get<decltype(Ns)>();
  int D=lBox.size();

  if (j.find("initialConditionGenerator") != j.end() )
    {
      auto & confGenJ = j["initialConditionGenerator"];
      if (confGenJ.find("lBox") != confGenJ.end() )
	{
	  lBoxInitialCondition=confGenJ["lBox"].get<decltype(lBoxInitialCondition)>() ;
	}
      else
	{
	  lBoxInitialCondition=lBox;
	}
	
    }
  else
    {
      lBoxInitialCondition = lBox;
    }
  if ( ( D != getDimensions() )or (lBoxInitialCondition.size() != getDimensions() ) )
    {
      throw invalidInput("Input file implies dimensionality different from " + std::to_string(getDimensions() ) );
      
    }
  
  
  geometryPBC geo( lBox[0], lBox[1], lBox[2]);
  
  states_t states;
  
  for (int i=0;i<Ns.size();i++)
    {
      state_t particleData(Ns[i] , D);
      particleData.setRandom();

      for (int i=0;i<getN(particleData);i++)
	{
	  for (int d=0;d<getDimensions();d++)
	    {
	      particleData(i,d)*=lBoxInitialCondition[d]*0.5;
	    }
	}
      //particleData*=lBox[0]/2.;
      states.push_back(particleData);
    }
  
  getFactory().registerJastrow< gaussianJastrow >();
  getFactory().registerJastrow< jastrowSquareWell >();
  getFactory().registerOrbital<sinOrbital>();
  getFactory().registerOrbital<planeWave>();
  getFactory().registerObservable<pairCorrelation>();
  getFactory().registerObservable<centerOfMassSquared>();
  getFactory().registerObservable<trimerhyperRadius>();
  
  auto waves = getFactory().createWavefunctions( j["wavefunctions"],geo);
  
  
  productWavefunction psi(waves);
  
  /* Potentials
     First register implemented concrete potential classes and then create the local potential as sum of individual potentials
*/

  getFactory().registerPotential<harmonicPotential>();
  getFactory().registerPotential<squareWellPotential2b>();
  
  auto potentials = getFactory().createPotentials(j["potentials"],geo);
  
  sumPotentials pot(potentials);
  
  auto eO=new  energy(&pot);
  auto efO= new forceEnergy(&pot);
  
  auto m = new realScalarEstimator("energy",eO);
  auto m2= new realScalarEstimator("forceEnergy",efO);
  
  std::string method = j["method"];

  real_t timeStep = j["timeStep"];
  size_t stepsPerBlock = j["stepsPerBlock"];
  size_t nBlocks = j["nBlocks"];  
  int seed= j["seed"];
  size_t correlationSteps=j["correlationSteps"];
  
  
  std::vector<states_t> configurations;

  if (j.find("configurations") != j.end())
    {
      configurations=readStatesFromDirectory(j["configurations"]);
      if (
	  (configurations.size() == 0)  and
	  ( (j.find("initialConfigurations") != j.end()) )

	  )
	{
	  configurations=readStatesFromDirectory(j["initialConfigurations"]);
	}
      
      if(  (pTools::rank() == 0 ) and ( configurations.size() == 0) )
	{
	  std::cout << ansiColor("yellow") << "WARNING: Failed to load configurations from directory. Falling back to randomly generated configurations." << ansiColor("default") << std::endl;
	}
      
    }

  auto ests = getFactory().createEstimators(j["measurements"]);
  auto storers = getFactory().createStorers(j["measurements"]);
  
  
  if (method == "vmc")
    {
      vmcDriver vmcO(&psi,timeStep);
      vmcO.getStepsPerBlock()=stepsPerBlock;
      vmcO.getCorrelationSteps()=correlationSteps;
      

      vmcO.getEstimators().push_back(m);
      
      if ( find(j["measurements"],"kind","forceEnergy") != j["measurements"].end() )
	{
	  vmcO.getEstimators().push_back(m2);
	}
      
      for (auto & est : ests)
	{
	  vmcO.getEstimators().push_back(est); 
	}
      
      
      vmcO.getRandomGenerator().seed(seed + pTools::rank() );
      states_t * initialConfiguration = &states;
      if (configurations.size() > 0 )
	{
	  initialConfiguration = &(configurations[0]);
	}

      if (! check_n_particles(*initialConfiguration,Ns) )
	{
	  throw invalidInput("Initial configuration does not math the numper of particles defined in the input file");
	}
      
      vmcO.run(*initialConfiguration,nBlocks); 
    }
  
  else if ( method == "dmc" or method == "svmc")
    {
      size_t nW=j["walkers"];
      
      dmcDriver dmcO(&psi,&pot,timeStep,nW);
      dmcO.getStepsPerBlock()=stepsPerBlock;
      dmcO.getCorrelationSteps()=correlationSteps;
      
      std::vector<states_t> dmcStates(configurations);      
      
      if ( dmcStates.size() == 0 )
	{
	  for(int i=0;i<nW/pTools::nProcesses();i++)
	    {
	      dmcStates.push_back(states);
	    }
      
	}
      
      dmcO.getRandomGenerator().seed(seed + pTools::rank() );
      
      if (method == "svmc")
	{
	  dmcO.disableBranching();
	  
	  if ( find(j["measurements"],"kind","forceEnergy") != j["measurements"].end() )
	{
	  dmcO.getEstimators().push_back(m2);
	}
		
	}

      
      for (auto & st : storers)
	{
	  dmcO.getEstimators().push_back(st);
	}
      
      
      for (auto & est : ests)
	{
	  dmcO.getEstimators().push_back(est);
	}
      

      dmcO.run(dmcStates,nBlocks);
    }

  pTools::finalize();
  
}
