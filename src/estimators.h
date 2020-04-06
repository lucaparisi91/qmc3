#ifndef ESTIMATORS_H
#define ESTIMATORS_H


#include <iostream>
#include "traits.h"
#include "accumulators.h"
#include "observables.h"
#include <iostream>
#include <fstream>     

class productWavefunction;
class walker;

class estimatorBase
{
public:
	using wave_t = productWavefunction;
	using walker_t=walker;
	
	estimatorBase(std::string label);
	virtual void accumulate(walker_t &w, wave_t & psi )=0;
	const std::string & label() const{ return _label;};
	virtual void write(std::ostream & stream)=0;
	virtual void clear() = 0;
	virtual void dump();

	virtual ~ estimatorBase();
private:
	std::string _label;
	std::fstream f;
	std::string filename;
};

template<class observable_t>
class estimator : public estimatorBase
{
public:
	using accumulator_t = typename observable_t::accumulator_t ;
	estimator(std::string label,observable_t  * ob_) : 
		estimatorBase::estimatorBase(label),ob(ob_) {}
	virtual void accumulate(walker_t & w,wave_t & psi) {ob->accumulate(w,psi,acc);}
	virtual void write(std::ostream & stream);
	virtual void clear(){acc.clear();}

private:
	accumulator_t acc;
	observable_t *ob;
};

class realScalarEstimator : public estimator<realScalarObservable>
{
public:
	using estimator<realScalarObservable>::estimator;
};



#endif