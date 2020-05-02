#include "estimatorsFactory.h"


realHistogramEstimator* createEstimatorFromOb(realHistogramObservable * ob,const json_t & j)
{
  return new realHistogramEstimator(ob,j);
}
