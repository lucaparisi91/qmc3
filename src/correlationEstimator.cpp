#include "correlationEstimator.h"
#include "walkers.h"
#include "tools.h"

realScalarStorer::realScalarStorer(std::string label_,realScalarObservable * ob_, int recordSteps_) : ob(ob_),recordSteps(recordSteps_),storer::storer(label_) {}

realScalarStorer::realScalarStorer(realScalarObservable * ob_,const json_t & j )  : realScalarStorer(j["label"].get<std::string>(),ob_,j["recordSteps"].get<int>() + 1 ) {};

void realScalarStorer::reserve(realScalarStorer::walker_t & w)
{
  w.getStorageScalarCorrelators()[getLabel()].resize(recordSteps);
  w.getTimeIndex()[getLabel()]=0.;
}

void realScalarStorer::reset(walker_t & w)
{
  w.getTimeIndex().at(getLabel() )=0.;
}

void realScalarStorer::store( walker_t & w, wavefunction_t & psi )
  {
    int & i = w.getTimeIndex().at(getLabel() );
 
    w.getStorageScalarCorrelators().at(getLabel() )(i)=(*ob)(w,psi);
    i=(i+1)% recordSteps;
  }

void realScalarForwardWalkingEstimator::accumulate(walker_t & w,wavefunction_t & psi)
  {
    const auto & i = w.getTimeIndex().at(targetLabel);
    
    auto & data =w.getStorageScalarCorrelators().at(targetLabel);
    auto recordSteps = data.size();
    assert(forwardWalkingSteps < recordSteps);
    auto j =  wrapIndex(i - forwardWalkingSteps - 1 ,recordSteps ) ;
    const auto & v = data( j);
    
    

    getAccumulator()+=v;
    
  }

realScalarForwardWalkingEstimator::realScalarForwardWalkingEstimator(std::string label,std::string targetLabel_, int forwardWalkingSteps_) : estimator<realScalarAccumulator_t>(label),targetLabel(targetLabel_),forwardWalkingSteps(forwardWalkingSteps_)
{
  
  std::ifstream of;
  of.open(getFileName());
  if (is_empty(of) & pTools::rank() == 0 )
    {
      auto & f=getFileDescriptor();
      f << getLabel() << std::endl;
    }
  of.close();
  
}

realScalarForwardWalkingEstimator::realScalarForwardWalkingEstimator( const json_t & j) : realScalarForwardWalkingEstimator( j["label"] ,  j["targetLabel"].get<std::string>(), j["forwardWalkingSteps"].get<int>() )
{
  
}
