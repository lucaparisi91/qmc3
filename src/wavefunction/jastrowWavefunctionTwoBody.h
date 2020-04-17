#include "wavefunction/wavefunction.h"

template<class jastrow_t>
class jastrowTwoBodyWavefunctionIndistinguishible : public wavefunction
{
public:
  
  jastrowTwoBodyWavefunctionIndistinguishible(jastrow_t J_,const geometry_t  &geo_, int setA_=0,int setB_=0) : setA(setA_), setB(setB_),J(J_),wavefunction::wavefunction(geo_)
  {
    if (setA == setB) throw invalidInput("setA == setB in distinguishable two body jastrow");
  }
  
 virtual real_t operator()(const walker_t & walker) 
	{		
	  auto & dis = walker.getTableDistances().distances(setA,setB);
	  
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
    auto & stateA = walker.getStates()[setA];
    auto & stateB = walker.getStates()[setB];

    auto & gradientA = walker.getGradients()[setA];
    auto & gradientB = walker.getGradients()[setB];
    
    auto & laplacian = walker.getLaplacianLog();
    auto & distances = walker.getTableDistances().distances(setA,setB);
    auto & differences = walker.getTableDistances().differences(setA,setB);
    auto & waveValue = walker.getLogWave();
    
    
    const int NA = getN(stateA);
    const int NB=getN(stateB);    
    constexpr int D = getDimensions();
    real_t  tmp,tmp1,tmp2;
    
    int k=0;
    for (int i=0;i<NA;i++)
      {
	for (int j=0;j<NB;j++)
	  {
	    auto d = distances(k);
	    
	    J.evaluateDerivatives(d,tmp,tmp1,tmp2);
	    
	    laplacian+=tmp2 + (D-1)*tmp1/d;
	    waveValue+=tmp;
	
	    for(int id=0;id<D;id++)
	      {
		gradientA(i,id)+=differences(k,id)/d * tmp1;
		gradientB(j,id)-=differences(k,id)/d * tmp1;
	      }
				
	    k++;
	  }
		
      
      }
  }

  
    
private:
    
  int setA;
  int setB;
  jastrow_t J;

};
