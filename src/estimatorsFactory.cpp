#include "estimatorsFactory.h"


realHistogramEstimator* createEstimatorFromOb(realHistogramObservable * ob,const json_t & j)
{
  return new realHistogramEstimator(ob,j);
}

realScalarEstimator* createEstimatorFromOb(realScalarObservable * ob,const json_t & j)
{
  return new realScalarEstimator(ob,j);
}


template<>
struct observableTraits<realScalarObservable >
{
  using storer_t = realScalarStorer;
  using estimator_t = realScalarEstimator;
  
};
