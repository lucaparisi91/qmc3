#include <random>
#include "traits.h"

class metropolis
{
public:
  metropolis(std::ranlux24 * rand) : _rand(rand),uniformDis(0.,1.) {}
  bool acceptLog(real_t ratioLog)
  {
    if ( ratioLog > 0. )
      {
	     return true;
      }
    else
      {
	     if ( ratioLog > std::log(uniformDis(*_rand)) )
	     {
	    return true;
	     }
	     else
	     {
	     return false;
	     }
      }
  }
  
private:
  std::ranlux24 *_rand;
  std::uniform_real_distribution<real_t> uniformDis;
  
};