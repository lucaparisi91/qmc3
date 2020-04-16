#include "orbitals.h"
#include "slaters.h"


double getLog(double a,int & sign){if (a < 0 ) {sign*=-1;};return std::log(std::abs(a));}
std::complex<double> getLog(std::complex<double> a,int & sign){return std::log(a*(1.*sign));sign=1;}

void tableSlaters::add(int setA,orbitalSetBase * orbitalSet)
{
  auto index=orbitalSets.size();
  orbitalSets.push_back(orbitalSet);
  indices1b[setA]=index;
  slaterMatrices.resize(slaterMatrices.size()+1);
  slaterMatricesInverse.resize(slaterMatricesInverse.size()+1);
  logDeterminants.resize(logDeterminants.size()+1);
  signs.resize(signs.size()+1);  
}

void tableSlaters::update(const tableSlaters::states_t & states)
{
  for ( const auto  element : indices1b ) // updates single set slaters
	{
	  auto & matrix=slaterMatrices[element.second];
	  auto & matrixInverse=slaterMatricesInverse[element.second];
	  auto & logDeterminant= logDeterminants[element.second];
	  auto & sign = signs[element.second];
	  
	  orbitalSets[element.second]->storeEvaluate(states[element.first],matrix);
	  /* compute inverse */
	  lud.compute(matrix);
	  matrixInverse=lud.inverse();
	  
	  /* Compute determinant*/
	  logDeterminant=0;
	  auto & luMatrix=lud.matrixLU();
	  sign=lud.permutationP().determinant();
	  for(auto i=0;i<luMatrix.rows();i++)
	    {
	      logDeterminant+=getLog(luMatrix(i,i),sign );
	    }
	  
	  logDeterminants[element.second]=logDeterminant;
	  

	}
}
