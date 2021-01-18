#include "toolsPimc.h"
#include <cmath>

using Real = double;

class poschlTellerPotential3D
{
    public:

    poschlTellerPotential3D(const json_t & j) : 
    poschlTellerPotential3D::poschlTellerPotential3D(j["V0"].get<Real>(),j["R"].get<Real>()  ) {}

    poschlTellerPotential3D(Real V0_,Real R_) : V0(V0_),R(R_) {}
    Real operator()(Real x,Real y, Real z) { Real r = std::sqrt(x*x + y*y + z*z); return  -V0/std::pow(cosh(r/R),2);}

    Real gradX(Real x,Real y , Real z) {Real r = std::sqrt(x*x + y*y + z*z); return x/r *radialDerivative(r) ;   }

    Real gradY(Real x,Real y , Real z) {Real r = std::sqrt(x*x + y*y + z*z); return y/r *radialDerivative(r) ;   }

    Real gradZ(Real x,Real y , Real z) {Real r = std::sqrt(x*x + y*y + z*z); return z/r *radialDerivative(r) ;   }

    

    private:

    Real radialDerivative(Real r) { Real x=r/R; return 2*V0/R*sinh(x)/(std::pow(cosh(x),3)) ;}


    Real V0;
    Real R;

};