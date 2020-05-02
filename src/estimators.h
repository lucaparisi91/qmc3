#ifndef ESTIMATORS_H
#define ESTIMATORS_H


#include <iostream>
#include "traits.h"
#include "accumulators.h"
#include "observables.h"
#include <iostream>
#include <fstream>     
#include <memory>


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
  virtual std::string & getLabel() {return _label;}
  virtual std::string & getFileName() {return filename;}

  virtual void accumulateMPI(int root)=0;
  virtual ~ estimatorBase();
  virtual std::vector<int> sets() const = 0;
  
  virtual std::fstream & getFileDescriptor() {return f;}
  
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

  virtual void accumulateMPI(int root) {acc.accumulateMPI(root);}
  auto & getAccumulator() {return acc;}
  
  virtual std::vector<int> sets() const {return ob->sets();};
  
private:
  accumulator_t acc;
  std::unique_ptr<observable_t> ob;
};

class realScalarEstimator : public estimator<realScalarObservable>
{
public:
  
  
  realScalarEstimator(std::string label_,realScalarObservable * ob_);
  using estimator<realScalarObservable>::estimator;
};

class realHistogramEstimator : public estimator<realHistogramObservable>
{
public:


  
  realHistogramEstimator(std::string label,realHistogramObservable * ob_,size_t size,real_t minx,real_t maxx);

  realHistogramEstimator(realHistogramObservable * ob_,const json_t & j);
  
  virtual void write(std::ostream & stream) override;
private:
  std::vector<real_t> x;
};

#endif
