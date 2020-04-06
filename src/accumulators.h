#ifndef ACCUMULATORS_H
#define ACCUMULATORS_H

#include "traits.h"


template<class T>
class scalarAccumulator
{
public:
	using value_t = T;
	scalarAccumulator() : sum(0.),n(0) {}

	void operator+=(value_t e){sum+=e;n+=1;};

	value_t average() const {return sum/n;}
	void clear(){sum=0.;n=0.;}
private:
	value_t sum;
	size_t n;
};

template<class T>
class vectorAccumulator
{
public:
	using vec_t = Eigen::Tensor<T,1> ;
	using value_t = T ;
	vectorAccumulator(){};

	vectorAccumulator(size_t size) : _sums(size),ns(size) {_sums.setConstant(0);ns.setConstant(0);}

	void accumulate(value_t a, size_t i) {ns[i]+=a;ns[i]+=1;}

	vec_t average() {return sums/ns;}

	const vec_t & sums() const {return _sums;}

	size_t size(){return  _sums.size();}
	
private:
	T _sums;
	Eigen::Tensor< real_t, 1> ns;
};


using realScalarAccumulator_t = scalarAccumulator<real_t>;
using realVectorAccumulator_t = vectorAccumulator<real_t>;

#endif