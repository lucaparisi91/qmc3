#include "driver.h"
#include "wavefunction/productWavefunction.h"
#include "estimators.h"
#include "walkers.h"
#include "initializer.h"
#include <iostream>
#include "tools.h"
namespace vmc
{

 vmcDriver::vmcDriver(wavefunction_t * wave_,real_t sigma_) :
 wave(wave_),
 sigma(sigma_),distribution(0.,1.) ,
 metropolisObj(&randGenerator),
 iBlock(0),nAccepted(0)
 {

 }

void update(walker_t & w,wavefunction_t & psi)
	{
		w.getTableDistances().update(w.getStates());
		w.getLogWave()=psi(w.getStates());
	};


void vmcDriver::accumulate()
{
	for( auto & est : _estimators)
	{
		est->accumulate(current_walker,*wave);
	}
}

void vmcDriver::run(states_t & states,size_t nBlocks)
{
	// initialize walkers
	initializer::initialize(current_walker,states,*wave);
	initializer::initialize(old_walker,states,*wave);
	initializer::initialize(tmp_walker,states,*wave);


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


void vmcDriver::step()
{
	std::swap(old_walker,current_walker);
	// update particle positions
	 auto & current_states = current_walker.getStates();
	const auto & old_states = old_walker.getStates();

	for (int i=0;i<current_states.size();i++)
	{
		auto & current_state = current_states[i];
		const auto & old_state = old_states[i];

		int N= current_state.dimensions()[0];
		int D= current_state.dimensions()[1];
		for(int i=0;i<N;i++)
		for (int d=0;d<D;d++)
		{
			current_state(i,d)=old_state(i,d) + distribution(randGenerator)*sqrt(sigma); 
		};

	}
	// update distances and evaluates the wavefunction
	update(current_walker,*wave);
	// accept or reject the walker
	bool accept = metropolisObj.acceptLog(2* (current_walker.getLogWave() - old_walker.getLogWave()));
	if (accept)
	{
		nAccepted+=1;
	}
	else
	{
		current_walker=old_walker;
	}
}

void vmcDriver::out()
{
	std::cout << ansiColor("green") << "Block: "<< ansiColor("default")<<iBlock<<std::endl;
	std::cout << "Acc. Ratio: " << nAccepted * 1./ _stepsPerBlock << std::endl;
	// erase
	auto energyEst = estimators()[0];
	std::cout << ansiColor("cyan") << "Energy: " << ansiColor("default");
	energyEst->write(std::cout);
	std::cout << std::endl;

	std::cout << std::endl;	


	nAccepted=0;
	 auto & ests = estimators();
	 ests.dump();
	 ests.clear();



}

}