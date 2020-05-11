#include "hyperRadius.h"
#include "walkers.h"
trimerhyperRadius::trimerhyperRadius(int setA_,int setB_,int setC_) : setA(setA_),setB(setB_),setC(setC_)
{
  
}

trimerhyperRadius::trimerhyperRadius(const json_t & j) : trimerhyperRadius(j["sets"][0],j["sets"][1],j["sets"][2])
{
  
}

real_t trimerhyperRadius::operator()(walker_t & w,wavefunction_t & psi)
{
  auto & norms1 =w.getTableDistances().distances(setA,setC);
  auto & norms2 =w.getTableDistances().distances(setB,setC);

  auto & norms3 = w.getTableDistances().distances(setA,setB);
  real_t radius=0;

  auto N1 = norms1.size();
  auto N2 = norms2.size();
  auto N3 = norms2.size();
  
  for (int i=0;i<N1;i++ )
    for (int j=0;j<N2;j++ )
      for (int k=0;k<N3;k++ )
	{
	  radius+=std::sqrt( norms1(i)*norms1(i) + norms2(j)*norms2(j) + norms3(k)*norms3(k) );
	}
  
  radius=radius/(N1*N2*N3);
  
  return radius;
  
}
