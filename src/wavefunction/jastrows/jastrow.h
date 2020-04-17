#ifndef JASTROW_H
#define JASTROW_H

#include <vector>
#include "parameters.h"
#include <iostream>
#include "qmcExceptions.h"
#include <nlohmann/json.hpp>


template<class T>
class jastrow
{
	/*
	CRT base class for the jastrow functions. Add default funcionalities to derived classes. 
	Derived class must posses at least a d0(x),d1(x)and d2(x) twhiche compute the zero,first and second derivative of the jastrow function
	*/
public:

	void registerParameter(int sourceParameter,parameter param)
	{
		sourceParameters.push_back(sourceParameter);
		parameters.push_back(param);

		if (sourceParameter >= static_cast<T*>(this)->nParameters())
		{
			throw invalidInput("Source Parameter larger than the number of parameters of the jastrow");
		};

	}
	 void addGradientParameters(real_t x,std::vector<real_t> & gradientParameter)
	 {
	 	for (int i=0;i<sourceParameters.size();i++)
	 		{
				static_cast<T*>(this)->addGradientParameter(x, sourceParameters[i], parameters[i], gradientParameter);
	 		}
	 }

	 void evaluateDerivatives(real_t x, real_t & d0_,real_t & d1_, real_t & d2_) {d0_=static_cast<T*>(this)->d0(x);d1_=static_cast<T*>(this)->d1(x);d2_=static_cast<T*>(this)->d2(x);}

	 void addGradientParameter(real_t x,int sourceParameter,parameter & param, std::vector<real_t> & gradientParameter) {throw missingImplementation("Jastrow does not seem to support any parameter gradient.");} ;

	 int nParameters() {return 0;}
  

protected:
	jastrow(){}; // disalloes the instantation of base class. Only concrete jastrows shoulf derive from this class

	std::vector<int> sourceParameters;
	std::vector<parameter> parameters; 
};

class gaussianJastrow : public jastrow<gaussianJastrow>
{
	// J(x) = -alpha* x^2 
public:
  int nParameters() {return 1;}; // number of variational parameters supported
  
  gaussianJastrow(real_t alpha_) : alpha(alpha_){};
  gaussianJastrow(const nlohmann::json  & j)
  {
    alpha=j["alpha"];
  }

  
  static std::string name() {return "gaussian"; }

  
  
  real_t d0(real_t x) {return -alpha*x*x;}
  real_t d1(real_t x) {return -2.*alpha*x;}
  real_t d2(real_t x) {return -2*alpha;}


	void addGradientParameter(real_t x,int sourceParameter,parameter & param, std::vector<real_t> & gradientParameter) 
	{
		if (sourceParameter == 0)
			{
				*(param.begin(gradientParameter) )-=x*x; 
			}
	};

	
	private:
	real_t alpha;
};

#endif
