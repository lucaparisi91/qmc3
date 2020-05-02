#include "estimators.h"
#include <fstream>
#include "tools.h"
#include "ptools.h"


template<class observable_t>
void estimator<observable_t>::write(std::ostream & stream)
	{
		stream << acc.average() ;
	};
template class estimator<realScalarObservable>;

estimatorBase::estimatorBase(std::string label) : _label(label)
{
	filename= label+".dat";
	
	f.open(filename,std::fstream::out | std::fstream::app);
}

estimatorBase::~estimatorBase()
{
	f.close();
}

void estimatorBase::dump()
{
	write(f);
	f << std::endl;
}

realScalarEstimator::realScalarEstimator(std::string label_,realScalarObservable * ob_) : estimator<realScalarObservable>::estimator(label_,ob_)
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

realHistogramEstimator::realHistogramEstimator(std::string label,realHistogramObservable * ob_,size_t size,real_t minx,real_t maxx) : estimator<realHistogramObservable>::estimator(label,ob_)
{
  getAccumulator()=accumulator_t(size,minx,maxx);
  
  std::ifstream of;
  of.open(getFileName());
  if (is_empty(of) & pTools::rank() == 0 )
    {
      auto & f=getFileDescriptor();
      f << "x " << getLabel() << std::endl;
    }
  
  of.close();

  x.resize(size);
  for (int i=0;i<x.size();i++)
    {
      x[i]=getAccumulator().stepSize() * i;
    }
  
}

realHistogramEstimator::realHistogramEstimator(realHistogramObservable * ob_,const json_t & j): realHistogramEstimator(j["label"].get<std::string>(),ob_,j["bins"].get<size_t>(), j["minx"].get<real_t>(),j["maxx"].get<real_t>()   )
{
  
}


void realHistogramEstimator::write(std::ostream & stream)
{
  auto out=getAccumulator().average();
  
  for (int i=0;i<x.size();i++)
    {
      stream << x[i] << " " <<  out[i] << std::endl ;
    }
};
