#include "walkers.h"
#include "dmcDriver.h"
#include "wavefunction/productWavefunction.h"
#include "initializer.h"
#include "tools.h"
#include "moves/dmcMoves.h"
#include "estimators.h"

bool noMetropolisPolicy::accept(
		mover & move,
		const dmcWalker & new_walker,
		const dmcWalker & old_walker,
		wavefunction_t & wavefunction,
		randomGenerator_t &generator)
{
	return true;
}

real_t noMetropolisPolicy::getAcceptanceRatio() const
{
  return 1.;
}


bool metropolisPolicy::accept(
		mover & move,
		const dmcWalker & new_walker,
		const dmcWalker & old_walker,
		wavefunction_t & wavefunction,
		randomGenerator_t &generator)
{
  real_t transitionRatio = move.transitionProbabilityRatio(new_walker,old_walker);

  return metropolisSampler.acceptLog(2*(new_walker.getLogWave() - old_walker.getLogWave() )  + transitionRatio ,generator);
  
}


real_t metropolisPolicy::getAcceptanceRatio() const
{
  return metropolisSampler.getAcceptanceRatio();
}

void metropolisPolicy::clear() 
{
  return metropolisSampler.clear();
}



 void dmcDriver::step()
{
	std::swap(current_walkers,old_walkers);
	current_walkers.resize(old_walkers.size());
	
	for (int i=0;i<old_walkers.size();i++)
	{
	  
		dmcMover->move(current_walkers[i],old_walkers[i],getRandomGenerator());
		updateForceGradientEnergy(current_walkers[i], getWavefunction(),energyOb);
		
		bool accepted=accepter->accept(*dmcMover,current_walkers[i],old_walkers[i],getWavefunction() ,getRandomGenerator());
		
		if (!accepted)
		{
		  current_walkers[i]=old_walkers[i];
		}

	}
}

void dmcDriver::run( const std::vector<states_t> &states , size_t nBlocks )
{
  current_walkers.resize(states.size());
  old_walkers.resize(states.size());
  
  for (int i=0;i<current_walkers.size();i++)
    {
      initializer::initialize(current_walkers[i],states[i],getWavefunction());
      initializer::initialize(old_walkers[i],states[i],getWavefunction());
    }
  
  driver::run(nBlocks);
}

void dmcDriver::out()
{
  	std::cout << ansiColor("green") << "Block: "<< ansiColor("default")<<getCurrentBlock()<<std::endl;
	
	auto & energyEst  = getEstimators()[0];
	std::cout << "Acc. Ratio: " << accepter->getAcceptanceRatio() << std::endl;
	std::cout << ansiColor("cyan") << "Energy: " << ansiColor("default");
	std::cout << std::scientific;
	energyEst->write(std::cout);
	std::cout << std::endl<<std::defaultfloat;
	
	auto & ests = getEstimators();
	ests.dump();
	ests.clear();

}

dmcDriver::dmcDriver(dmcDriver::wavefunction_t * wave, potential * pot,real_t timeStep) :
  driver::driver(wave), energyOb(pot),
  dmcMover( new driftDiffusionFirstOrder(timeStep)),
  energyEst(new realScalarEstimator("energy",&energyAccFromWalker) ),
  accepter(new metropolisPolicy)
{
  getEstimators().push_back( energyEst.get() );
}


void dmcDriver::accumulate()
{
	auto & wave=getWavefunction();

	for( auto & est : getEstimators())
	{
	  for (auto & current_walker : current_walkers)
		est->accumulate(current_walker,wave);
	}
}
