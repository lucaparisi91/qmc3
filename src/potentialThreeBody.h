#ifndef POTENTIAL_THREE_BODY_H
#define POTENTIAL_THREE_BODY_H


#include "potential.h"
#include "tools3B.h"

template<class functor_t>
class potentialThreeBodyUnDis :  public potential
{
public:
  using potential::operator();
  
  potentialThreeBodyUnDis( const functor_t & V_,int setA_,const geometry_t & geo) : potential(geo),V(V_),setA(setA_){};
  
  potentialThreeBodyUnDis(const json_t & j, const geometry_t & geo) :
    potentialThreeBodyUnDis(functor_t(j["functor"]),j["sets"][0],geo)
  {
    
  }

  
  virtual real_t operator()(const walker_t & w)
  {
    const auto & dis=w.getTableDistances().distances(setA,setA);
    real_t sum=0;

    const int N = getN(w.getStates()[setA]);

    
    LOOP3B( N,
	   auto rji = dis(ji);	
	   auto rkj = dis(kj);	
	   auto rki = dis(ki);
	    
	   auto R= sqrt(rji*rji + rki*rki + rkj*rkj);
	    
	   sum+=V(R);
	   
	   )

    return sum;
    
  }
  
  static std::string name() {return "potentialThreeBodyUnDis/" + functor_t::name(); }
  
  virtual std::vector<int> sets() const {return {setA,setA,setA};} ;
  
  
private:
  
  functor_t V;
  int setA;
  
};


#endif
