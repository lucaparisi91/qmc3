#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include "traits.h"
#include <unsupported/Eigen/CXX11/Tensor>
#include "geometry.h"
#include "wavefunction/jastrows/jastrow.h"
#include "qmcExceptions.h"


struct wavefunctionComponentCommands;



class wavefunction
{
	/*
	Represents the logarithm of the wavefunction. Supports wavefunctions acting on at most two different sets
	*/
public:
	using geometry_t = geometry;
	using state_t = Eigen::Tensor<real_t, 2>;
	using grad_t = state_t;
	using states_t=std::vector<state_t>;
	using grads_t = std::vector<grad_t>;

	wavefunction(const geometry_t & geo_ );
	wavefunction(const geometry_t & geo_ , int setA_);


	virtual real_t operator()(const state_t & state ) {throw missingImplementation("Evaluation on a single species.");};
	virtual real_t operator()(const state_t & state1,const state_t & state2 ) {throw missingImplementation("Evaluation on dual species.");};

	const auto & getGeometry() {return *geo;}

	void evaluateDerivatives(const states_t & state, grads_t & gradient , real_t & wavevalue, real_t & laplacian);


	virtual real_t operator()(const states_t & state );

	virtual void evaluateDerivatives(const state_t & state, grad_t & gradient , real_t & wavevalue, real_t & laplacian) {throw missingImplementation("Derivatives evaluation on single species.");}; // evaluates all derivatives and the value of the wavefunction in one go for efficiency 
	virtual void evaluateDerivatives(const state_t & state1, const state_t & state2,const grad_t & gradient1 , const grad_t & gradient2, const real_t & wavevalue, const real_t & laplacian) {throw missingImplementation("Derivatives evaluation on two species.");}; // evaluates all derivatives and the value of the wavefunction in one go for efficiency . Gradients are added to the input vectors

	//virtual void evaluateDerivatives(const tableDistances & state,grad_t & gradient,real_t & wavevalue,real_t & laplacian) {throw missingImplementation("Derivatives evaluation on single species from distance table.");}; // evaluates all derivatives and the value of the wavefunction in one go for efficiency 

	~wavefunction();
private:
	const geometry_t * geo;
	wavefunctionComponentCommands * commands;

};


struct wavefunctionComponentCommands
{
	using states_t= wavefunction::states_t;
	using grads_t=wavefunction::grads_t;

	virtual real_t operator()(const states_t & states){throw missingImplementation("Evaluation on a vector of components");return 0;}; 
	virtual void evaluateDerivatives(const states_t & state, grads_t & gradient , real_t & wavevalue, real_t & laplacian){throw missingImplementation("Evaluation on a vector of components");}
};

struct wavefunctionSingleComponentCommands : wavefunctionComponentCommands
{
	wavefunctionSingleComponentCommands(wavefunction * w_,int setA_) : setA(setA_),w(w_){};
	using states_t= wavefunction::states_t;
	virtual real_t operator()(const states_t & states){return (*w)(states[0]);}
	virtual void evaluateDerivatives(const states_t & state, grads_t & gradient , real_t & waveValue, real_t & laplacian){(*w).evaluateDerivatives(state[0], gradient[0],waveValue,laplacian);}
private:
	wavefunction * w;
	int setA;

};



template<class jastrow_t>
class jastrowOneBodyWavefunction :  public wavefunction
{
public:
	using diff_t = Eigen::Tensor<real_t,2>;
	using distances_t= Eigen::Tensor<real_t,1>;

	jastrowOneBodyWavefunction(jastrow_t J_,const geometry_t  &geo_) : J(J_),wavefunction::wavefunction(geo_) {}
	jastrowOneBodyWavefunction(jastrow_t J_,const geometry_t  &geo_, int setA) : J(J_),wavefunction::wavefunction(geo_,setA) {}
	virtual real_t operator()(const state_t & state) 
	{
		differences = this->getGeometry().differencesOneBody(state,{0,0,0});
		distances = norm(differences);
		exit(0);
		int N=state.dimensions()[0];
		real_t sum=0;

		for(int i=0;i<N;i++)
		{
			sum+=J.d0(distances(i));
		}
		return sum;
	};

	virtual void evaluateDerivatives(const state_t & state, grad_t & gradient , real_t & waveValue, real_t & laplacian)
	{
		real_t tmp,tmp1,tmp2;
		differences = this->getGeometry().differencesOneBody(state,{0,0,0});
		distances = norm(differences);
		int N=state.dimensions()[0];
		int D=state.dimensions()[1];

		laplacian=0;
		for (int i=0;i<N;i++)
			{
				auto d = distances(i);
				J.evaluateDerivatives(d,tmp,tmp1,tmp2);
				laplacian+=tmp2 + (D-1)*tmp1/d;
		 		waveValue+=tmp;
				for(int id=0;id<D;id++)
					{
	  					gradient(i,id)+=differences(i,d)*(tmp1/d);
					}
			}
      
	}

private:
	jastrow_t J;
	diff_t differences; // temporary variable
	distances_t distances; // temporary variable
	std::array<real_t, 3> center;
};

#endif