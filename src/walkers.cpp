#include "walkers.h"
#include "energy.h"
#include "wavefunction/productWavefunction.h"

void updateForceGradientLaplacian(walker & w,productWavefunction & psi)
{
	/* Update forces ,laplacian and wavefunction value*/
  w.getTableDistances().update(w.getStates());
  w.getTableSlaters().update(w.getStates());
  
  psi.evaluateDerivatives(w);
};

void updateForceGradientEnergy(dmcWalker & w,productWavefunction & psi, energy & energyOb)
{
  w.getTableSlaters().update(w.getStates());
  w.getTableDistances().update(w.getStates());
  w.getEnergy()=energyOb(w,psi);
};

void dmcWalker::createMPIDataType()
{
  /*
    Create a data type which send particle configurations, force gradients , laplacian and energy. Does NOT send slater and distances information ( which in standarard dmc are not required for the evolution)
*/
  const int N = getStates().size();
  const int M= 2;
  const int NScalars= 3;
  const int Ntot = N*M + NScalars;
  MPI_Aint  offsets[Ntot] ;
  int blockCounts[Ntot];
  MPI_Datatype dtypes[Ntot];
  int k=0;
  for (int i=0;i<N;i++)
    {
      MPI_Get_address(getStates()[i].data(), &offsets[k]);
      blockCounts[k]=getStates()[i].size();
      dtypes[k]=MPI_DOUBLE;
      k++;
    }
  
  for (int i=0;i<N;i++)
    {
      MPI_Get_address(getGradients()[i].data(), &offsets[k]);
      blockCounts[k]=getStates()[i].size();
      dtypes[k]=MPI_DOUBLE;
    }
  k++;
  MPI_Get_address(&getLogWave(), &offsets[k]);
  blockCounts[k]=1;
  dtypes[k]=MPI_DOUBLE;
  k++;
  MPI_Get_address(&getLaplacianLog(), &offsets[k]);
  blockCounts[k]=1;
  dtypes[k]=MPI_DOUBLE;

  k++;
  MPI_Get_address(&getEnergy(), &offsets[k]);
  blockCounts[k]=1;
  dtypes[k]=MPI_DOUBLE;

  

  
  

  
  MPI_Type_struct(Ntot, blockCounts, offsets, dtypes, &getMPIDatatype() );
  MPI_Type_commit(&getMPIDatatype());
  
  //MPI_Get_address(getStates().data(), &offsets[0]);
  // MPI_Get_address(&getLogWave() , &offsets[1] );
  // MPI_Get_address(&getGradients() , &offsets[2] );
  // MPI_Get_address(&getLaplacianLog() , &offsets[3] );
  // MPI_Get_address(&getEnergy() , &offsets[4] );

  

  
  
}
