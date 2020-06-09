#include "jastrow.h"
#include <cmath>


class jastrowHardSphereGauss : public jastrow<jastrowHardSphereGauss>
{
public:
  
  jastrowHardSphereGauss(real_t a_,real_t Rm_, real_t D_) ;

  jastrowHardSphereGauss(const json_t & j);
    
  inline real_t d0(const real_t & x) const {return x<=Rm ? scattering_d0(x) : gaussian_d0(x) ;}
  inline real_t d1(const real_t & x) const {return x<=Rm ? scattering_d1(x) : gaussian_d1(x) ;}
  inline real_t d2(const real_t & x) const {return x<=Rm ? scattering_d2(x) : gaussian_d2(x) ;}
  
  inline real_t gaussian_d0(real_t x) const {return logC  - alpha*(x-d)*(x-d)  ;}
  inline real_t gaussian_d1(real_t x) const {return -2* alpha*(x-d)  ;}
  inline real_t gaussian_d2(real_t x) const {return -2* alpha  ;}
  
  inline real_t scattering_d0(real_t x) const {return log(1 - a/x);}
  inline real_t scattering_d1(real_t x) const {return a/(x*(x-a));}
  inline real_t scattering_d2(real_t x) const  {return -a/std::pow(x*(x-a),2) * (2*x - a); }
  
  static std::string name() {return "hardSphereGauss";};
  
  
  
  private:    
  real_t Rm;
  real_t d;
  real_t a;
  real_t C;
  real_t alpha;
  real_t logC;
    
};
