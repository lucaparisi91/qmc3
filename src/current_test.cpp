#include "parameters.h"
#include "wavefunction/jastrows/jastrow.h"
#include <iostream>
#include <unsupported/Eigen/CXX11/Tensor>
#include "geometry.h"
#include "wavefunction/productWavefunction.h"
#include "tableDistances.h"
#include "potential.h"
#include "initializer.h"
#include "walkers.h"
#include "energy.h"
#include "estimators.h"
#include "driver.h"
#include "dmcDriver.h"
#include "branching.h"
#include <nlohmann/json.hpp>
#include "centerOfMassSquared.h"
#include "wavefunction/jastrowWavefunctionOneBody.h"
#include "wavefunction/jastrowWavefunctionTwoBody.h"

#include "wavefunction/slaterDeterminant.h"

#include "orbitals.h"

#include "wavefunction/jastrows/jastrowSpline.h"


auto getEnergies(const   pTools::walkerDistribution::walkers_t & walkers)
{
  /* Gather all energies on process 0 for this walker distribution */
  std::vector<real_t> energies;
  for (int i=0;i<walkers.size();i++)
    {
      energies.push_back(walkers[i].getEnergy() );
    }

  std::vector<int> populations;
  populations.resize(pTools::nProcesses() , 0 );
  int n = energies.size();
  MPI_Gather( & n,1 , MPI_INT,
	      populations.data(), 1, MPI_INT,
	      0, MPI_COMM_WORLD);

  
  std::vector<real_t> totalEnergies;
  int totN = 0;
  for (int i=0;i<populations.size();i++)
    {
      totN+=populations[i];
    }
  std::vector<int> offsets;
  offsets.resize(populations.size(),0);
  
  std::partial_sum(populations.begin(),populations.end() - 1,offsets.begin() + 1 );
  totalEnergies.resize(totN);
  
  MPI_Gatherv(
	      energies.data(),
	      energies.size(),
	      MPI_DOUBLE,
	      totalEnergies.data(),
	      populations.data(),
	      offsets.data(),
	      MPI_DOUBLE,
	      0,
	      MPI_COMM_WORLD
	      );
  return totalEnergies;
  
}



int main(int argc, char** argv)
{
  pTools::init(argc,argv);


  Eigen::ArrayXd coefficients(102);

  
  coefficients << -6.80202700e-05,  3.40101350e-05, -6.80202700e-05, -3.74111485e-04,
        -8.84263511e-04, -1.59847635e-03, -2.51674999e-03, -3.63908445e-03,
        -4.96547971e-03, -6.49593579e-03, -8.23045267e-03, -1.01690304e-02,
        -1.23116689e-02, -1.46583682e-02, -1.72091283e-02, -1.99639493e-02,
        -2.29228310e-02, -2.60857736e-02, -2.94527769e-02, -3.30238411e-02,
        -3.67989661e-02, -4.07781519e-02, -4.49613985e-02, -4.93487059e-02,
        -5.39400741e-02, -5.87355032e-02, -6.37349930e-02, -6.89385437e-02,
        -7.43461552e-02, -7.99578274e-02, -8.57735605e-02, -9.17933544e-02,
        -9.80172091e-02, -1.04445125e-01, -1.11077101e-01, -1.17913138e-01,
        -1.24953236e-01, -1.32197395e-01, -1.39645614e-01, -1.47297895e-01,
        -1.55154236e-01, -1.63214638e-01, -1.71479101e-01, -1.79947624e-01,
        -1.88620209e-01, -1.97496854e-01, -2.06577560e-01, -2.15862327e-01,
        -2.25351155e-01, -2.35044043e-01, -2.44940992e-01, -2.55042003e-01,
        -2.65347073e-01, -2.75856205e-01, -2.86569398e-01, -2.97486651e-01,
        -3.08607965e-01, -3.19933340e-01, -3.31462776e-01, -3.43196272e-01,
        -3.55133830e-01, -3.67275448e-01, -3.79621127e-01, -3.92170867e-01,
        -4.04924668e-01, -4.17882529e-01, -4.31044451e-01, -4.44410434e-01,
        -4.57980478e-01, -4.71754583e-01, -4.85732748e-01, -4.99914975e-01,
        -5.14301262e-01, -5.28891610e-01, -5.43686018e-01, -5.58684488e-01,
        -5.73887018e-01, -5.89293609e-01, -6.04904261e-01, -6.20718974e-01,
        -6.36737748e-01, -6.52960582e-01, -6.69387477e-01, -6.86018433e-01,
        -7.02853450e-01, -7.19892528e-01, -7.37135666e-01, -7.54582866e-01,
        -7.72234126e-01, -7.90089447e-01, -8.08148828e-01, -8.26412271e-01,
        -8.44879774e-01, -8.63551338e-01, -8.82426963e-01, -9.01506649e-01,
        -9.20790396e-01, -9.40278203e-01, -9.59970071e-01, -9.79866000e-01,
    -9.99965990e-01, -1.02027004e+00;

  
  real_t stepSize=0.010101010101010102;
  
  jastrowSpline j(coefficients,stepSize,-1,-2);
  
  real_t x=0.5;
  real_t d0,d1,d2;

  real_t dx=1e-4;
  
  j.evaluateDerivatives(x,d0,d1,d2);
  
  for (real_t x=0.1;x<1;x+=dx)
    {
      j.evaluateDerivatives(x,d0,d1,d2);

      std::cout << x << " " << d0 << " " << " " <<  d1 << " " << d2 << std::endl;    
    }
  
  pTools::finalize();

}
