#ifndef JASTROWWAVEFUNCTIONTHREEBODY_H
#define JASTROWWAVEFUNCTIONTHREEBODY_H

#include "wavefunction/wavefunction.h"

#include "tools3B.h"



template<class jastrow_t>
class jastrowThreeBodyWavefunctionUnDistinguishable : public wavefunction
{
public:
  
  jastrowThreeBodyWavefunctionUnDistinguishable(jastrow_t J_,const geometry_t  &geo_, int setA_=0) : setA(setA_),J(J_),wavefunction::wavefunction(geo_)
  {
    
  }

  
  jastrowThreeBodyWavefunctionUnDistinguishable(const json_t & j,const geometry_t & geo ) : jastrowThreeBodyWavefunctionUnDistinguishable( jastrow_t(j["jastrow"]), geo,j["sets"][0]  ) {}

  virtual real_t operator()(const walker_t & walker)
  {
    auto & dis = walker.getTableDistances().distances(setA,setA);  
  
    real_t sum=0;
    int N = getN( walker.getStates()[setA] );

    // performs a 3b loop on all particles
    
    LOOP3B( N ,
	    auto rji = dis(ji);	
	    auto rkj = dis(kj);	
	    auto rki = dis(ki);
	    
	    auto R= sqrt(rji*rji + rki*rki + rkj*rkj);
	    sum+=J.d0(R);
	    ) ;
  
    return sum;
    
  }

  
  virtual std::vector<int> sets() const {return {setA,setA} ;}
  
  virtual void accumulateDerivatives( walker_t & walker ) override
  {
    auto & state = walker.getStates()[setA];

    auto & gradient = walker.getGradients()[setA];
    
    auto & laplacian = walker.getLaplacianLog();
    auto & dis = walker.getTableDistances().distances(setA,setA);
    auto & diff = walker.getTableDistances().differences(setA,setA);
    auto & waveValue = walker.getLogWave();
    
    
    const int N = getN(state);
    constexpr int D = getDimensions();
    
    real_t  d0,d1,d2;
    
    
    LOOP3B(N,

	   
	   auto rji = dis(ji);	
	   auto rkj = dis(kj);	
	   auto rki = dis(ki);
	   
	   auto R= std::sqrt(rji*rji + rki*rki + rkj*rkj);
	   
	   J.evaluateDerivatives(R,d0,d1,d2);
	   
	   
	   for (int id=0;id < D;id++)
	     {    
	       gradient(i,id)-=  d1*( diff(ji,id) + diff(ki,id))/R;
	       gradient(j,id)+=  d1*(   diff(ji,id) - diff(kj,id) )/R;
	       gradient(k,id)+=  d1*(   diff(ki,id) + diff(kj,id) )/R;
	     }
	   
	   laplacian+= d2 * 3 ;
	   laplacian += 3 * d1 / R *  (2*D-1);
	   
	   waveValue+=d0;
	   
	   
	   

	   )
  }
  
  static std::string name()   {return "jastrow3bUnDis/" + jastrow_t::name();}
  
  virtual std::string print() const override {
    
    return J.print(0,getGeometry().getLBox(0)/2. , 10000)
      ;}  
  
private:
  int setA ;
  jastrow_t J;
  
};



#endif
