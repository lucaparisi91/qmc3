#include "metropolis.h"
#include "traits.h"
#include "walkers.h"
#include "estimatorCollection.h"


class productWavefunction;

class estimatorBase;
class walker;

namespace vmc
{
	using wavefunction_t=productWavefunction;
	using walker_t = walker;
void update(walker_t & w,wavefunction_t & psi);

class vmcDriver
{
public:
	using estimators_t = estimatorCollection;
	vmcDriver(wavefunction_t * wave,real_t sigma_);
	size_t & stepsPerBlock() {return _stepsPerBlock;}
	void run(states_t & states,size_t nSamples);
	void step(); // perform a mc step
	void accumulate(); // accumulate measurements
	void out(); // output summery of the blocks

	auto & estimators() {return _estimators;}

	auto & currentWalker() {return current_walker;}
	auto & oldWalker(){return old_walker;}
private:

	walker_t current_walker;
	walker_t old_walker;
	walker_t tmp_walker;

	metropolis metropolisObj;
	estimators_t _estimators;
	wavefunction_t *wave;
	real_t sigma;
	real_t nAccepted;
	std::ranlux24 randGenerator;
	std::normal_distribution<real_t> distribution;
	size_t _stepsPerBlock;
	size_t iBlock;

};


}