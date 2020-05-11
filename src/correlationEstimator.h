#ifndef CORRELATION_ESTIMATOR_H
#define CORRELATION_ESTIMATOR_H

#include "estimators.h"
#include "accumulators.h"
#include "qmcExceptions.h"


class unImplementedStorer;
template<class T>
struct estimatorTraits
{
  using storer_t = unImplementedStorer;
};

class realScalarStorer;

template <>
struct estimatorTraits<realScalarAccumulator_t>
{
  using storer_t = realScalarStorer;
};

class storer
{
public:
  using walker_t = walker;
  
  storer(std::string label_) : label(label_) {}
  virtual void reserve(walker_t & w) = 0;
  virtual void store( walker_t & w,wavefunction_t & psi )=0;
  const auto & getLabel() {return label;}
  virtual std::vector<int> sets() const {return {} ; };
private:  
  std::string label;
};

class unImplementedStorer : public storer
{
public:
  
  template<class T > unImplementedStorer(T * , const json_t & j) : storer("invalid"){throw missingImplementation("No storer class  was implemented.");}
   virtual void reserve(walker_t & w) {throw missingImplementation("No reserve was implemented.");};
  virtual void store( walker_t & w,wavefunction_t & psi ) {throw missingImplementation("No store was implemented.");};
  
};


class realScalarStorer : public storer
{
public:
  using observable_t=realScalarObservable;
  
  realScalarStorer(std::string label_,realScalarObservable * ob_, int recordSteps_) ;
  realScalarStorer(realScalarObservable * ob_,const json_t & j ) ;
  virtual void reset( walker_t & w );
  
  virtual void reserve(walker_t & w);
  
  virtual void store( walker_t & w, wavefunction_t & psi );

  virtual std::vector<int> sets() const override {return ob->sets();}
private:
  std::unique_ptr<observable_t> ob;
  int recordSteps;
};

class realScalarForwardWalkingEstimator : public estimator<realScalarAccumulator_t>
{
public:
  realScalarForwardWalkingEstimator(std::string label,std::string targetLabel_, int forwardWalkingSteps_) ;
  
  realScalarForwardWalkingEstimator(const json_t & j);
  
  virtual void accumulate(walker_t & w,wavefunction_t & psi) override;

  virtual void write(std::ostream & w) override;
private:

  std::string targetLabel;
  int forwardWalkingSteps;
  
};


#endif
