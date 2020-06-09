#include "walkers.h"
#include "jastrowWavefunctionThreeBody.h"



hardSphereConstraintThreeBodyUnDis::hardSphereConstraintThreeBodyUnDis(int setA_,real_t R0) : setA(setA_),V(R0) {}

hardSphereConstraintThreeBodyUnDis::hardSphereConstraintThreeBodyUnDis(const json_t & j) :
  hardSphereConstraintThreeBodyUnDis::hardSphereConstraintThreeBodyUnDis(j["sets"][0].get<real_t>(), j["hardRadius"].get<real_t>() )
{
  
}


bool hardSphereConstraintThreeBodyUnDis::operator()(const walker_t & w)
  {

    const auto & dis=w.getTableDistances().distances(setA,setA);
    
    const int N = getN(w.getStates()[setA]);
    
    LOOP3B( N,
	   auto rji = dis(ji);	
	   auto rkj = dis(kj);	
	   auto rki = dis(ki);
	    
	   auto R= sqrt(rji*rji + rki*rki + rkj*rkj);
	    
	    if (R <= V.R0() )
	      {
		return false;
	      }
	    )
      
    return true;  
  }

