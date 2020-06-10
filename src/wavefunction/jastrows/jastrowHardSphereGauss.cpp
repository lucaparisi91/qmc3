#include "jastrowHardSphereGauss.h"


jastrowHardSphereGauss::jastrowHardSphereGauss(real_t a_,real_t Rm_,real_t D_) : a(a_),Rm(Rm_),d(D_)
{
  alpha=a/(Rm*(Rm-a)*2*(d-Rm) );
  C=(1-a/Rm)*exp(alpha*(Rm-d)*(Rm-d) );
  logC=log(C);
  
};


jastrowHardSphereGauss::jastrowHardSphereGauss(const json_t & j) : jastrowHardSphereGauss(j["a"],j["Rm"],j["D"] )
{
  
};

jastrowHardSphere::jastrowHardSphere(real_t a_) : a(a_)
{
  
};

jastrowHardSphere::jastrowHardSphere(const json_t & j) : jastrowHardSphere(j["a"].get<real_t>() )
{
  
};

