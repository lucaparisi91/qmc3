#include "accumulators.h"
class productWavefunction;
class walker;

template<class T>
class observable
{
public:
	observable() {}
	using accumulator_t = T;
	using value_t = typename T::value_t ;
	using wavefunction_t = productWavefunction;
	using walker_t = walker;

	virtual void accumulate(walker_t & walker, wavefunction_t & wavefunction,  accumulator_t & acc)=0;	
};

class realScalarObservable : public observable<realScalarAccumulator_t> 
{
public:
	virtual void accumulate(walker_t & w, wavefunction_t & wavefunction,  accumulator_t & acc) override
	{
		acc+=(*this)(w,wavefunction);
	}

	virtual real_t operator()(walker_t & walker,wavefunction_t & wave)=0;
};
