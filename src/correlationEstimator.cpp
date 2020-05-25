#include "correlationEstimator.h"
#include "walkers.h"
#include "tools.h"

realScalarStorer::realScalarStorer(std::string label_,realScalarObservable * ob_, int recordSteps_) : ob(ob_),recordSteps(recordSteps_),storer::storer(label_) {}

realScalarStorer::realScalarStorer(realScalarObservable * ob_,const json_t & j )  : realScalarStorer(j["label"].get<std::string>(),ob_,j["recordSteps"].get<int>()  ) {};

void realScalarStorer::reserve(realScalarStorer::walker_t & w)
{
  w.getStorageScalarCorrelators()[getLabel()].resize(recordSteps);
  w.getTimeIndex()[getLabel()]=0.;
  w.getFillingStatus()[getLabel()]=true;
}

void realScalarStorer::reset(walker_t & w)
{
  w.getTimeIndex().at(getLabel() )=0.;
  w.getFillingStatus().at(getLabel())=true;
}

void realScalarStorer::store( walker_t & w, wavefunction_t & psi )
  {
    int & i = w.getTimeIndex().at(getLabel() );
    
    w.getStorageScalarCorrelators().at(getLabel() )(i)=(*ob)(w,psi);
    i=(i+1)% recordSteps;
    if (i==0)
      {
	w.getFillingStatus().at(getLabel())=false;
      }
  }

void realScalarForwardWalkingEstimator::accumulate(walker_t & w,wavefunction_t & psi)
  {
    const auto & i = w.getTimeIndex().at(targetLabel);
    
    auto & data =w.getStorageScalarCorrelators().at(targetLabel);
    auto recordSteps = data.size();
    assert(forwardWalkingSteps < recordSteps);
    auto j =  wrapIndex(i - forwardWalkingSteps - 1 ,recordSteps ) ;
    const auto & v = data( j);
    
    if ( ! w.getFillingStatus().at(targetLabel) )
      {
	getAccumulator()+=v;
      }
    
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

void realScalarForwardWalkingEstimator::write(std::ostream & stream)
{
    if ( getAccumulator().getWeight() > 0 )
    {
      estimator<realScalarAccumulator_t>::write(stream);
    }

}



realHistogramStorer::realHistogramStorer(std::string label_,realHistogramObservable * ob_,size_t size,real_t minx,real_t maxx,int recordSteps_) : ob(ob_),recordSteps(recordSteps_),storer::storer(label_) ,tmpAcc(size,minx,maxx)
{
  
}

void realHistogramStorer::store(realHistogramStorer::walker_t & w,wavefunction_t & psi)
{
  int & i = w.getTimeIndex().at(getLabel() );
  tmpAcc.clear();
  ob->accumulate(w,psi,tmpAcc);
  
  auto & data=w.getStorageScalarCorrelators().at(getLabel() );
  
  for(int j=0;j<tmpAcc.size();j++)
    {
      data(i *tmpAcc.size() + j)= tmpAcc.sums()(j)/tmpAcc.weight();
    }
  
  i=(i+1)% recordSteps;
  if (i==0)
    {
      w.getFillingStatus().at(getLabel())=false;
    }
}


void realHistogramStorer::reserve(walker_t & w)
{
  w.getStorageScalarCorrelators()[getLabel()].resize(recordSteps*tmpAcc.size() );
  w.getTimeIndex()[getLabel()]=0.;
  w.getFillingStatus()[getLabel()]=true;
  
}

void realHistogramStorer::reset(walker_t & w)
{
  w.getTimeIndex().at(getLabel() )=0.;
  w.getFillingStatus().at(getLabel())=true;
}

realHistogramStorer::realHistogramStorer(realHistogramObservable * ob_,const json_t & j ) : realHistogramStorer(j["label"],ob_,j["bins"],j["minx"],j["maxx"],j["recordSteps"])
{
  
}

void realHistogramForwardWalkingEstimator::accumulate(walker_t & w,wavefunction_t & psi)
  {
    const auto & i = w.getTimeIndex().at(targetLabel);

    auto & acc=getAccumulator();
    
    
    auto & data =w.getStorageScalarCorrelators().at(targetLabel);
    auto recordSteps = data.size()/getAccumulator().size();
    assert(forwardWalkingSteps*acc.size() < data.size());
    
    auto j =  wrapIndex(i - forwardWalkingSteps - 1 ,recordSteps ) ;
    
    auto & sums = acc.sums();
    if ( ! w.getFillingStatus().at(targetLabel) )
      {
	for(int ii=0;ii<acc.size();ii++)
	  {
	    sums(ii)+=data(j*acc.size() + ii );
	  }
	acc.weight()+=1;
      }
    
  }



realHistogramForwardWalkingEstimator::realHistogramForwardWalkingEstimator(std::string label,std::string targetLabel_,size_t size,real_t minx,real_t maxx,int forwardWalkingSteps_) :  estimator<realHistogramAccumulator_t>(label),targetLabel(targetLabel_),forwardWalkingSteps(forwardWalkingSteps_)
{
  getAccumulator().resize(size,minx,maxx);


  
  std::ifstream of;
  of.open(getFileName());
  if (is_empty(of) & pTools::rank() == 0 )
    {
      auto & f=getFileDescriptor();
      f << "x" <<" "<<getLabel() << std::endl;
    }
  
  of.close();
}


void realHistogramForwardWalkingEstimator::write(std::ostream & stream)
{
    if ( getAccumulator().weight() > 0 )
    {
      estimator<accumulator_t>::write(stream);
    }

}

realHistogramForwardWalkingEstimator::realHistogramForwardWalkingEstimator(const json_t & j) : realHistogramForwardWalkingEstimator(j["label"].get<std::string>(),j["targetLabel"].get<std::string>(),j["bins"].get<int>(),j["minx"].get<real_t>(),j["maxx"].get<real_t>(),j["forwardWalkingSteps"].get<int>() )
{
  
}
