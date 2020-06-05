#ifndef TRIMER_CORRELATION_H
#define TRIMER_CORRELATION_H

#include "observables.h"

class trimerCorrelation : public realHistogramObservable
{
public:
  trimerCorrelation(int setA,int setB,int setC);
  trimerCorrelation(const json_t & j);
  
  static std::string name() {return "trimerCorrelation";}
  
  virtual std::vector<int> sets() const override {return {setA,setB,setC}; }

  virtual void accumulate(walker_t & w,wavefunction_t & wavefunction,accumulator_t & acc);

private:

  void setNormalizationFactor(const walker_t & w , const wavefunction_t & psi ,const accumulator_t & acc) ;

  int setA;
  int setB;
  int setC;

  real_t _normalizationFactor;
  real_t max;
  real_t deltax;
  
};

#endif