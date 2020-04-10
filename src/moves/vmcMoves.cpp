#include "vmcMoves.h"
#include "walkers.h"
#include <iostream>

gaussianMover::gaussianMover(real_t sigma_ ) 
: sigma(sigma_),distribution(0.,1.)
{
  
}
void gaussianMover::move( 
	gaussianMover::walker_t & new_walker,
	const gaussianMover::walker_t & old_walker,
	gaussianMover::rand_t & randGenerator

	)
{
	
	auto & current_states = new_walker.getStates();
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
}
