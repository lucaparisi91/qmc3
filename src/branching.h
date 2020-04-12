#ifndef BRANCHING_H
#define BRANCHING_H
#include "traits.h"
class dmcWalker;

class branchingControl 
{
public:
  using walker_t = dmcWalker;
  branchingControl(real_t timeStep_,real_t meanWalkers_,real_t deltaWalkers_);

  void branch(std::vector<walker_t> & newWalkers,const std::vector<walker_t> & oldWalkers,randomGenerator_t & rand);
  int nDescendants(walker_t & new_walker, const walker_t & old_walker,randomGenerator_t & rand);

  void setEnergyShift(const std::vector<walker_t> & walkers);
  
private:
  size_t meanWalkers;
  size_t deltaWalkers;
  size_t maxWalkers;
  real_t energyShift; // E_T , can be offset to costrain the number of walkers
  real_t timeStep;
  std::uniform_real_distribution<real_t> uniformDis;
  std::vector<int> _nDescendants;
};

#endif
