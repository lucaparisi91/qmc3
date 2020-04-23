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
