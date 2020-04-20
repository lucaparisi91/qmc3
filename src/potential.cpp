#include "potential.h"
#include "geometry.h"
#include "tableDistances.h"
#include <iostream>
#include "tools.h"
#include "walkers.h"


potential::potential(const geometry_t & geo_ ) : geo(&geo_)
{

};

potential1b::potential1b(const geometry_t & geo_ , int setA_) : potential::potential(geo_),_setA(setA_)
{
	
};

harmonicPotential::harmonicPotential(const geometry_t & geo,real_t freq , int setA ) : potential1b(geo,setA) ,omega(freq) {}

real_t harmonicPotential::operator()(const walker_t & w)
{
  const auto & dis= w.getTableDistances().distances(setA() );

  real_t sum=0;
  
  for(int i=0;i<dis.size();i++)
    {
      sum+=0.5 * omega * omega * dis(i) * dis(i);
    }  
  return sum;
};

sumPotentials::sumPotentials(std::vector<potential*> potentials_) :
  _potentials(potentials_)
{
  
}

real_t sumPotentials::operator()(const walker_t & state )
{
  real_t sum=0;
  for (int i=0;i<size();i++)
    {
      sum+=(*(_potentials[i]))(state);
    }
  return sum;
}

harmonicPotential::harmonicPotential(const json_t & j, const geometry_t & geo) : harmonicPotential(geo, j["omega"], j["set"] ) {};

squareWellPotential2b::squareWellPotential2b(const geometry_t & geo,real_t V0_ , real_t R0_ , int setA_, int setB_ ) : R0(R0_),setA(setA_),setB(setB_),V0(V0_), potential(geo) {}

squareWellPotential2b::squareWellPotential2b(const json_t & j, const geometry_t & geo) : squareWellPotential2b(geo,j["V0"],j["R0"],j["sets"][0],j["sets"][1]) {}


real_t squareWellPotential2b::operator()(const walker_t & w)
{
  const auto & dis= w.getTableDistances().distances(setA,setB);
  
  auto nPairsInt=std::count_if(dis.data(),dis.data() + dis.size(), [&] (const real_t  &d) {return d<= R0;});
  
  return -nPairsInt*V0;
  
};
