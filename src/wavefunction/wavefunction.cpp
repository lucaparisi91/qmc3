#include "wavefunction.h"
#include "tableDistances.h"

wavefunction::wavefunction(const geometry_t & geo_ ) : geo(&geo_)
{
  
};

bool wavefunction::satisfyConstraints(const walker_t & state)
{
  for (auto & storedConstraint : constraints)
    {
      if ( not (*storedConstraint)(state) )
	{
	  return false;
	}
    }
  return true;
};


void wavefunction::addConstraint( std::unique_ptr<constraint> p)
{
  constraints.push_back( std::move(p) );
  
  
}
