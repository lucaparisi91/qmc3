#include "estimators.h"

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