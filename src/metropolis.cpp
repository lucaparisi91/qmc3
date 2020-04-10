#include "metropolis.h"

metropolis::metropolis(std::ranlux24 * rand) : _rand(rand),uniformDis(0.,1.),n(0),nAccepted(0.) {}

bool metropolis::acceptLog(real_t ratioLog)
  {
    bool result;
    if ( ratioLog > 0. )
      {
	     result=true;
      }
    else
      {
	     if ( ratioLog > std::log(uniformDis(*_rand)) )
	     {
	    result=true;
	     }
	     else
	     {
	     result=false;
	     }
      }
    if (result)
      {
	nAccepted+=1;
      }
    n+=1;
    return result;
  }


real_t metropolis::getAcceptanceRatio()
{
  return nAccepted*1./n;
}

void metropolis::clear()
{
  nAccepted=0;
  n=0;
}
