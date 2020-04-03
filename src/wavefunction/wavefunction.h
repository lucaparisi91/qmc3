#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include "traits.h"
#include <unsupported/Eigen/CXX11/Tensor>
#include "geometry.h"
#include "wavefunction/jastrows/jastrow.h"
#include "qmcExceptions.h"


struct wavefunctionComponentCommands;
class tableDistances;



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


	virtual real_t operator()(const states_t & state ) = 0;
	virtual real_t operator()(const states_t & state,tableDistances & tab ) = 0;


	const auto & getGeometry() const {return *geo;}


	virtual void evaluateDerivatives(const states_t & state, grads_t & gradient , real_t & wavevalue, real_t & laplacian)=0;

	virtual void evaluateDerivatives(const states_t & state, grads_t & gradient , real_t & wavevalue, real_t & laplacian,const tableDistances & tab)=0;

	virtual std::vector<int> sets() const = 0 ;



private:
	const geometry_t * geo;

};


class wavefunction1b : public wavefunction
{
public:
	wavefunction1b(const geometry_t & geo, int setA_); 

	virtual real_t operator()(const states_t & states ){return (*this)(states[_setA]);}
	virtual real_t operator()(const states_t & state,tableDistances & tab ) ;


	virtual void evaluateDerivatives(const states_t & states, grads_t & gradient , real_t & wavevalue, real_t & laplacian)
	{
		evaluateDerivatives(states[_setA], gradient[_setA], wavevalue, laplacian);
	}

	virtual void evaluateDerivatives(const states_t & states, grads_t & gradient , real_t & wavevalue, real_t & laplacian,const tableDistances & tab);



	virtual real_t operator()(const state_t & state);


	virtual void evaluateDerivatives(const state_t & states, grad_t & gradient , real_t & wavevalue, real_t & laplacian);

	virtual void evaluateDerivatives(const state_t & states, grad_t & gradient , real_t & wavevalue, real_t & laplacian, const difference_t  & diff,const distance_t & distance)=0;

	virtual real_t operator()(const state_t & state,const distance_t & dis) { throw missingImplementation("Evaluation of single state form distances");}

	const int & setA() const {return _setA;}

	std::vector<int>  sets() const {return {_setA};}

private:
	int _setA;
	distance_t distances;
	difference_t differences;
};

template<class jastrow_t>
class jastrowOneBodyWavefunction :  public wavefunction1b
{
public:
	using diff_t = Eigen::Tensor<real_t,2>;
	using distances_t= Eigen::Tensor<real_t,1>;

	using wavefunction1b::operator();
	using wavefunction1b::evaluateDerivatives;


	jastrowOneBodyWavefunction(jastrow_t J_,const geometry_t  &geo_, int setA=0) : J(J_),wavefunction1b::wavefunction1b(geo_,setA) {}

	virtual real_t operator()(const state_t & state,const distance_t & dis) override 
	{
		
		
		int N=state.dimensions()[0];
		real_t sum=0;

		for(int i=0;i<N;i++)
		{
			sum+=J.d0(dis(i));
		}
		return sum;
	};

	virtual void evaluateDerivatives(const state_t & state, grad_t & gradient , real_t & waveValue, real_t & laplacian , const difference_t & differences,const distance_t & distances) override
	{
		real_t tmp=0,tmp1=0,tmp2=0;
		
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
	  					gradient(i,id)+=differences(i,id)/d * tmp1;
					}
				

			}
		
      
	}

private:
	jastrow_t J;
};

#endif