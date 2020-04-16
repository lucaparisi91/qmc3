#include "wavefunction.h"

template<class jastrow_t>
class jastrowOneBodyWavefunction :  public wavefunction
{
public:
  using diff_t = ::difference_t;
  using distances_t= ::difference_t;

  using wavefunction::accumulateDerivatives;
  
  jastrowOneBodyWavefunction(jastrow_t J_,const geometry_t  &geo_, int setA_=0) : setA(setA_), J(J_),wavefunction::wavefunction(geo_) {}
  
  virtual real_t operator()(const walker_t & walker) 
	{		
	  auto & dis = walker.getTableDistances().distances(setA);
	  
	  int N=size(dis);
	  real_t sum=0;

	  for(int i=0;i<N;i++)
	    {
	      sum+=J.d0(dis(i));
	    }
	  return sum;
	};
  
	virtual void accumulateDerivatives( walker_t & walker ) override
	{
	  auto & state = walker.getStates()[setA];
	  auto & gradient = walker.getGradients()[setA];
	  auto & laplacian = walker.getLaplacianLog();
	  auto & distances = walker.getTableDistances().distances(setA);
	  auto & differences = walker.getTableDistances().differences(setA);
	  auto & waveValue = walker.getLogWave();

	  
	  const int N = getN(state);
	  constexpr int D = getDimensions();
	  real_t  tmp,tmp1,tmp2;
	  
	 
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

  virtual std::vector<int> sets() const {return {setA};}

private:
  jastrow_t J;
  int setA;
};

