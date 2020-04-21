#ifndef PTOOLS_H
#define PTOOLS_H

#include <mpi.h>
#include <string>
#include <vector>
#include <stdint.h>
#include <limits.h>
#include "traits.h"


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

  
  int isend(state_t,int fromRank,int toRank);
  
  

  
  void determineLoadBalanceComunicationsAliasMethod( std::vector<int> & populations,std::vector<int> & permutations,std::vector<int> & sources,std::vector<std::vector<int> > & destinations,std::vector<int> & amounts)  ; // all operatans are modified

  

  template<class T>  class walkerContainer;

  class dmcWalker;

  
  
class walkerDistribution
{
public:
  using walkers_t = walkerContainer<dmcWalker>;
  
  walkerDistribution();
  void gatherPopulations(walkers_t & walkers); // allgather the populations  distribued across all processors
  
  void determineComm( const std::vector<int> & populations); // internally figure out the communications between the processors. Done redundantly on al processors
  
  void isend(const walkers_t & walkers );
  void ireceiv(const walkers_t & walkers);
  
  void wait(); // stall until walkers have been received

private:
  int _currentRank;
  int _nProcesses;
  std::vector<int> _sources;
  std::vector<int> _permutations;
  std::vector< std::vector<int> > _sendToRanks;
  std::vector<int> tmpPopulations;
  std::vector<int> nWalkersReceived;
  
};

};
#endif
