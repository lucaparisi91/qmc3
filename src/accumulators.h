#ifndef ACCUMULATORS_H
#define ACCUMULATORS_H

#include "traits.h"
#include "ptools.h"
#include <iostream>
template<class T>
class scalarAccumulator
{
public:
	using value_t = T;
	scalarAccumulator() : sum(0.),n(0) {}

	void operator+=(value_t e){sum+=e;n+=1;};

	value_t average() const {return sum/n;}
	void clear(){sum=0.;n=0.;}

  void accumulateMPI(int root) // sum partial sums on all processors into the root processor
  {
    
    sum=pTools::sum(sum,root);
    n=pTools::sum(n,root);
  }
private:
	value_t sum;
	size_t n;
};

template<class T>
class vectorAccumulator
{
public:
  using vec_t = Eigen::Matrix<T,Eigen::Dynamic,1> ;
  using value_t = T ;
  vectorAccumulator(){};

  vectorAccumulator(size_t size) : _sums(size),ns(size) {_sums.setConstant(0);ns.setConstant(0);}

  void accumulate(value_t a, size_t i) {ns[i]+=a;ns[i]+=1;}

  vec_t average() {return sums/ns;}

  const vec_t & sums() const {return _sums;}
  
  size_t size(){return  _sums.size();}

  
  
private:
  vec_t _sums;
  Eigen::VectorXd ns;
};


template<class T>
class histogramAccumulator
{
public:
  using vec_t = Eigen::Matrix<T,Eigen::Dynamic,1> ;
  using value_t = T ;
  histogramAccumulator(){};
  histogramAccumulator(size_t size,real_t min_, real_t max_) : _sums(size),_minx(min_),_maxx(max_),_weight(0) {_sums.setConstant(0);deltax=(_maxx-_minx)/size;deltaxInverse=1./deltax;}
  
  void accumulate(value_t a,value_t x) {int i=int((x-_minx)*deltaxInverse); _sums[i]+=a;}
  
  vec_t average() const {return sums()/weight();}
  
  const vec_t & sums() const {return _sums;}
  
  size_t size() const {return  _sums.size();}

  real_t maxx() const {return _maxx;}
  
  auto & weight() {return _weight;}
  const auto & weight() const {return _weight;}

  real_t stepSize() const {return deltax;}

  void clear(){_sums.setConstant(0.);_weight=0.;}

  void accumulateMPI(int root) // sum partial sums on all processors into the root processor
  {
    
    //_sums=pTools::sum(_sums,root);
    _weight=pTools::sum(_weight,root);
  }
  
private:
  vec_t _sums;
  T _weight;
  real_t _maxx;
  real_t _minx;
  real_t deltax;
  real_t deltaxInverse;
};


using realScalarAccumulator_t = scalarAccumulator<real_t>;
using realVectorAccumulator_t = vectorAccumulator<real_t>;
using realHistogramAccumulator_t = histogramAccumulator<real_t>;
#endif
