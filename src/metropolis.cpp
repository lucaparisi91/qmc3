#include "metropolis.h"

metropolis::metropolis() : n(0),nAccepted(0){}

bool metropolis::acceptLog(real_t ratioLog,randomGenerator_t & randGen)
  {
    bool result;
    if ( ratioLog > 0. )
      {
	     result=true;
      }
    else
      {
	if ( ratioLog > std::log(uniformDis(randGen)) )
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

real_t metropolis::getAcceptanceRatio() const
{
  return nAccepted*1./n;
}

void metropolis::clear()
{
  nAccepted=0;
  n=0;
}
