

class optimizationObservables : public realVectorObservable
{
public:
  
  optimizationObservables(std::vector<std::vector<mappedOptimizationParameter> > & parameters);
  
  virtual void accumulate(walker_t & w,wavefunction_t & wavefunction,accumulator_t & acc);
  
  optimizationObservables(const json_t & j);
  
  static std::string name() {return "optimizationObservables";}
  
private:
  
  std::vector<real_t> parameterGradient;
  
};
