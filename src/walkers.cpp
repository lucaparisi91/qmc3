#include "walkers.h"
#include "energy.h"
#include "wavefunction/productWavefunction.h"
#include "tools.h"
#include <sys/stat.h>

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








template<class T>
void walkerContainer<T>::push_back(const T &  w)
    {
      _size=_size +1;
    
      if (_size > capacity() )
	{
	  walkers.push_back(std::unique_ptr<T>() );
	  (*(walkers.end() -1 )).reset( new T(w));
	}
      else
	{
	  *(walkers[_size-1])=w;
	}
    }



template<class T>
void walkerContainer<T>::dump(int i)
  {
    std::ofstream f;
    int pId =pTools::rank();

    struct stat st = {0};

    if (stat(baseDir.c_str(), &st) == -1) {
      mkdir(baseDir.c_str(), 0700);
}
    
    f.open(baseDir + "/walkers-Rank" + std::to_string(pId) + ".json");

    json_t j = toJson();
    f << j;
    f.close();
    
  }

template<class T>
json_t walkerContainer<T>::toJson()
  {
    json_t j;
    j["configurations"]=json_t::array({});
    for (int i=0;i<walkers.size();i++)
      {
	j["configurations"].push_back(::toJson((*this)[i].getStates()) );
      }
    return j;
  }

template<class T>
void walkerContainer<T>::reserve(size_t size2,const T & w)
  {
    auto oldCap = capacity();
    if (size2 > capacity() ) walkers.resize(size2);
    for (int i=oldCap ; i < capacity() ;i++)
      {
	walkers[i].reset(new T(w));
      }
  }


template<class T>
void walkerContainer<T>::resize(size_t size2)
  {
    if (size2 > capacity()  )
      {
	
	resize(size2,T());
      }
    else
      {
	_size=size2;
      }
  }

template<class T>
void walkerContainer<T>::resize(size_t size2, const T & w)
  {
    reserve(size2,w);    
    _size=size2;
  }






#include "walkers.hpp"
