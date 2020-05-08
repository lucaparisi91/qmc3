#include "jastrowSquareWell.h"
#include <cmath>

void jastrowSquareWell::initCoefficients()
{
   K0=std::sqrt(V0);
   
   // as=R0*(1-tan(K0*R0)/(K0*R0));
   //std::cout << a << std::endl;
   
   halflBox=lBox/2.;
   //C=as/(Rm*Rm*alpha*(-exp(-alpha*Rm) + exp(-alpha*(lBox-Rm))));
   
   //B=1-as/Rm - C*(exp(-alpha*Rm) + exp(-alpha*(lBox-Rm)));
   
   A=std::sin(K0*R0)/(1-R0*aInverse);
   
   C=-A/(Rm*Rm*alpha* ( std::exp(-alpha*(lBox - Rm)) - std::exp(-alpha*Rm)  ));
   
   B=A*(1/Rm-aInverse) - C*( std::exp(-alpha*Rm) + std::exp(-alpha*(lBox-Rm)));
   
   // std::cout << "A " << A <<  std::endl;
   //std::cout << "K0 " << K0 <<  std::endl;
   // std::cout << "R0 " << R0 <<  std::endl;
   // std::cout << "B" << B << std::endl;
   // std::cout << "C" << C << std::endl;
    
   if (B + C*2*exp(-alpha*halflBox) < 0)
     {
       std::cout << "Jastrow cannot be negative" << std::endl;
       exit(1);
     }
   else
     {
       longDistanceConstant=std::log(B + C*2*std::exp(-alpha*halflBox));
     }
   
}


jastrowSquareWell::jastrowSquareWell(const json_t & j) :
  jastrowSquareWell(j["V0"],j["R0"],j["Rm"],j["alpha"],j["aInverse"],j["cut_off"].get<double>()*2){}
