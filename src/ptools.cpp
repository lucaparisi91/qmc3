#include "ptools.h"

#include <iostream>

namespace pTools
{

    int init(int argc,char** argv) {return MPI_Init(&argc,&argv);}
  
  int finalize(){return MPI_Finalize();}

  int nProcesses() {int res;  MPI_Comm_size(MPI_COMM_WORLD, &res); return res; }

  int rank() {int res;  MPI_Comm_rank(MPI_COMM_WORLD, &res); return res; }
  
  bool isMaster()
  {
    if (rank()==0)
      return true;
    else
      return false;
  }
  

  int broadcast(int * data, int count,int root)
  {
    return MPI_Bcast( data, count,MPI_INT,root,MPI_COMM_WORLD);
  }
  int broadcast(char * const data, int count,int root)
  {
    return MPI_Bcast( data, count,MPI_CHAR,root,MPI_COMM_WORLD);
  }

  
  int broadcast(int * data,int root)
  {
    return broadcast(data,1,root);
  }

  int broadcast(std::string * data,int root)
  {
    int size=data->size();
    broadcast(&size,root);
    data->resize(size);
    return broadcast(const_cast<char*>(data->data()),size,root); 
  }

  
  int sum(double * send,double * recv ,int count,int root)
  {
    return MPI_Reduce( send, recv,count, MPI_DOUBLE,MPI_SUM,root,MPI_COMM_WORLD);
  };

  double sum(const double & sum, int root )
  {
    double tmp=sum;
    int status = MPI_Reduce (&sum,&tmp,1,MPI_DOUBLE,MPI_SUM,root,MPI_COMM_WORLD);
    return tmp;
  }

  
  int sum(const int & sum, int root )
  {
    int tmp=sum;
    int status = MPI_Reduce (&sum,&tmp,1,MPI_INT,MPI_SUM,root,MPI_COMM_WORLD);
    return tmp;
  }

  size_t sum(const size_t & sum, int root )
  {
    size_t tmp=0;
    int status = MPI_Reduce (&sum,&tmp,1,MPI_SIZE_T,MPI_SUM,root,MPI_COMM_WORLD);
    return tmp;
  }
  
}
