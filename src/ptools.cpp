#include "ptools.h"
#include <iostream>
#include <numeric>
#include "walkers.h"

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

  void determineLoadBalanceComunicationsAliasMethod( std::vector<int> & populations,std::vector<int> & permutations,std::vector<int> & sources,std::vector<std::vector<int> > & destinations,std::vector<int> & amounts)
  {
  /* 
Sources: filled with the rank of the processor it must receive from
amounts: filled with number of walkers send by the source
populations : current population distribution
  */
  
  int size = std::accumulate(populations.begin(),populations.end(),0);
  permutations.resize(size);
  std::iota (std::begin(permutations), std::end(permutations), 0);
  // sort such that A[i] < k all on the left of A[j] > k
  int k= size/populations.size();
  // add phantom tasks
  int r=size%populations.size();
  
  int i=0;
  while( (i < populations.size())  and ( populations[i] < k) )
    {
      i+=1;
    }
  
  int j=i+1;
  
  while (j<populations.size() )
    {
      if (populations[j] < k )
	{
	  std::swap(populations[i],populations[j]);
	  std::swap(permutations[i],permutations[j]);
	  i+=1;
	}
      j+=1;
      
    }
  
  // alias algorithm to determine the amount of comunications
  
  j=0;
  sources.resize(populations.size(),0);
  amounts.resize(populations.size(),0);
  destinations.resize(populations.size(),{});  
  
  int I,J;
  auto balanced =  [&] ( int idx) {return  k + (1 ? idx < r : 0 );  } ;
  
  while ( i> j)
    {
      J = permutations[j];
      I = permutations[i];
      
      sources[J]=I;
      
      amounts[J]=balanced(j)- populations[j];
      populations[i]=populations[i] - amounts[J];
      destinations[I].push_back(J);
      
      if (populations[i] < balanced(i) )
	{
	  ++i;
	}
      j++;
    }
  }  


walkerDistribution::walkerDistribution()
{
  _currentRank=rank();
  _nProcesses=nProcesses();
  
}
  
void walkerDistribution::determineComm(const std::vector<int> & populations)
{
  tmpPopulations=populations;
  pTools::determineLoadBalanceComunicationsAliasMethod( tmpPopulations, _permutations, _sources, _sendToRanks, nWalkersReceived);
    }

  int isend(double * p,int count,int dest,int tag,MPI_Request &req)
  {
    return MPI_Isend( p ,  count,MPI_DOUBLE, dest, tag,
			  MPI_COMM_WORLD, &req);
  }

  
  int isend(state_t & state,int toRank,int tag,MPI_Request &req)
  {
    auto data_ptr = state.data();
    return isend(data_ptr,state.size(),toRank,tag,req);
  };

  int partialSend(dmcWalker & w,int destination,int tag)
  {
    w.createMPIDataType();
    return MPI_Send(MPI_BOTTOM, 1, w.getMPIDatatype(), destination, tag, MPI_COMM_WORLD);
  }
  
  int partialRecv(dmcWalker * w, int source,int tag)
  {
    w->createMPIDataType();
    MPI_Status status;
    
    return MPI_Recv(MPI_BOTTOM, 1, w->getMPIDatatype(), source, tag, MPI_COMM_WORLD, & status);

  }
};

