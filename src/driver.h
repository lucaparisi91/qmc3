#ifndef DRIVER_H
#define DRIVER_H

#include "metropolis.h"
#include "traits.h"
#include "walkers.h"
#include "estimatorCollection.h"
#include "moves/move.h"

class productWavefunction;

class estimatorBase;
class walker;

class driver
{
public:
    using estimators_t = estimatorCollection;
	using wavefunction_t = productWavefunction;
	driver(wavefunction_t * wave_);

	auto & getRandomGenerator() {return _rand;}

	size_t & getStepsPerBlock()  {return _stepsPerBlock;}

	size_t getCurrentBlock() const {return iBlock;}

	size_t getCurrentSubStep() const  {return iSubStep;}

	auto & getEstimators() {return _estimators;}

	virtual void run(size_t nBlocks);
	virtual void step()=0;
	virtual void accumulate()=0;
	virtual void out()=0;

	wavefunction_t & getWavefunction () {return (*_wave);}
	const wavefunction_t & getWavefunction() const {return (*_wave);}

private:
	wavefunction_t * _wave;
	std::ranlux24 _rand;
	size_t _stepsPerBlock;
	size_t iBlock;
	size_t iSubStep;
	estimators_t _estimators;

};




void update(walker & w,productWavefunction & psi);

class vmcDriver : public driver
{
public:
	using walker_t = walker;

	vmcDriver(wavefunction_t * wave,real_t sigma_);
	
	void run(states_t & states,size_t nSamples);
	void step(); // perform a mc step
	void out(); // output summery of the blocks

	virtual void accumulate() ; // accumulate measurements

	auto & currentWalker() {return current_walker;}
	auto & oldWalker(){return old_walker;}

private:
	walker_t current_walker;
	walker_t old_walker;
	walker_t tmp_walker;
	metropolis metropolisObj;
	std::unique_ptr<mover> vmcMove;
};


#endif
