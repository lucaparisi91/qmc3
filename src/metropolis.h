#ifndef METROPOLIS_H
#define METROPOLIS_H


#include <random>
#include "traits.h"

class metropolis
{
public:
  metropolis(std::ranlux24 * rand);
  bool acceptLog(real_t ratioLog);

  real_t getAcceptanceRatio();

  void clear();
private:
  std::ranlux24 *_rand;
  std::uniform_real_distribution<real_t> uniformDis;
  size_t n;
  size_t nAccepted;
};


#endif
