
#include "abstractFactory.h"
#include "estimators.h"
#include "estimatorCollection.h"


/*
Defines a wavefunction factory. Does not anything about jastrows. Id are strings formed by concateneting recursively 'kind' all kind items in the object using '/' as a separator.
*/


typedef estimatorBase* (*estimatorCreatorFunc) ( const json_t & j);

realHistogramEstimator* createEstimatorFromOb(realHistogramObservable * ob,const json_t & j);


template<class observable_t>
estimatorBase * createEstimator(const json_t & j )
{
  observable_t * observable = new observable_t(j);

  auto est= createEstimatorFromOb(observable,j);

  return est;
}


class estimatorFactory : public abstractFactory<estimatorBase,std::string, estimatorCreatorFunc>
{
public:
  using abstractFactory_t= abstractFactory<estimatorBase,std::string, estimatorCreatorFunc>;

  
  template<class ob_t >
  void registerObservable()
  {
    registerType( ob_t::name()  , & (createEstimator<ob_t> ) );
  }

  auto create(const json_t & j)
  {
    std::vector<estimatorBase*> estimators;

    for (auto & estJson : j )
      {
      
    

	std::string id = estJson["kind"];
	
	if ( (id != "forceEnergy") and (id !="energy") )
	  estimators.push_back( abstractFactory_t::create(id,estJson) );
      }
    
    return estimators;
  }
  
 
  
};
