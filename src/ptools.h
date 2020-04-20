#ifndef PTOOLS_H
#define PTOOLS_H


#include <mpi.h>
#include <string>

#include <stdint.h>
#include <limits.h>

#if SIZE_MAX == UCHAR_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
   #define MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
   #error "what is happening here?"
#endif



namespace pTools
{
  int init(int argc,char** argv);
  
  int finalize();

  int nProcesses() ;

  int rank() ;
  
  bool isMaster();  

  
  int broadcast(int * data, int count,int root);
  int broadcast(char * const data, int count,int root);  
  int broadcast(int * data,int root);  
  int broadcast(std::string * data,int root);

  int sum(double * ob,int count,int root);

  int sum(double * send,double * recv ,int count,int root);

  double sum(const double & sum,int root);
  size_t sum(const size_t & sum, int root );
  int sum(const int & sum, int root );

}

#endif
