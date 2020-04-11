#include "driver.h"
#include "wavefunction/productWavefunction.h"
#include "estimators.h"
#include "walkers.h"
#include "initializer.h"
#include <iostream>
#include "tools.h"
#include "moves/vmcMoves.h"


driver::driver(driver::wavefunction_t * wave) : _wave(wave), 
iBlock(0),iSubStep(0),_stepsPerBlock(0)
{

}

void driver::run(size_t nBlocks)
{

	for (iBlock=0;iBlock<nBlocks;iBlock++)
	{
		for (size_t i=0;i<_stepsPerBlock;i++)
		{
			step(); 
			accumulate();			
		};
		out();
	}
};




 vmcDriver::vmcDriver(vmcDriver::wavefunction_t * wave_,real_t sigma_) :
driver::driver(wave_),
 metropolisObj(),
 vmcMove(new gaussianMover(sigma_))
 {
 	
 }

void vmcDriver::run(states_t & states,size_t nBlocks)
{
	wavefunction_t & wave= getWavefunction();

	initializer::initialize(current_walker,states,getWavefunction());
	initializer::initialize(old_walker,states,wave);
	initializer::initialize(tmp_walker,states,wave);
	driver::run(nBlocks);

};

void update(walker & w,productWavefunction & psi)
	{
	  w.getTableDistances().update(w.getStates());
	  w.getLogWave()=psi(w.getStates());
	};


void vmcDriver::step()
{
	auto & wave=getWavefunction();

	std::swap(old_walker,current_walker);

	vmcMove->move(current_walker,old_walker,getRandomGenerator());

	// update distances and evaluates the wavefunction
	update(current_walker,wave);
	// accept or reject the walker
	bool accept = metropolisObj.acceptLog(2* (current_walker.getLogWave() - old_walker.getLogWave()),getRandomGenerator());
	if (!accept)
	{
		current_walker=old_walker;
	}
}

void vmcDriver::accumulate()
{
	auto & wave=getWavefunction();

	for( auto & est : getEstimators())
	{
		est->accumulate(current_walker,wave);
	}
}

void vmcDriver::out()
{
	std::cout << ansiColor("green") << "Block: "<< ansiColor("default")<<getCurrentBlock()<<std::endl;
	std::cout << "Acc. Ratio: " << metropolisObj.getAcceptanceRatio() << std::endl;
	// erase
	
	auto energyEst = getEstimators()[0];
	std::cout << ansiColor("cyan") << "Energy: " << ansiColor("default");
	energyEst->write(std::cout);
	std::cout << std::endl;
	
       
	auto & ests = getEstimators();
	ests.dump();
	ests.clear();
	metropolisObj.clear();
	 
}
