#include "estimatorsFactory.h"

realHistogramEstimator* createEstimatorFromOb(realHistogramObservable * ob,const json_t & j)
{
  return new realHistogramEstimator(ob,j);
}

realScalarEstimator* createEstimatorFromOb(realScalarObservable * ob,const json_t & j)
{
  return new realScalarEstimator(ob,j);
}

realVectorEstimator* createEstimatorFromOb(realVectorObservable * ob,const json_t & j)
{
  return new realVectorEstimator(ob,j);
}


template<>
struct observableTraits<realScalarObservable >
{
  using storer_t = realScalarStorer;
  using estimator_t = realScalarEstimator;
  
};


template<>
struct observableTraits<realHistogramObservable >
{
  using storer_t = realHistogramStorer;
  using estimator_t = realHistogramEstimator;
  
};


std::vector<estimatorBase*> estimatorFactory:: create(const json_t & j)
  {
    std::vector<estimatorBase*> estimators;
    
    for (auto & estJson : j )
      {

	std::string id = estJson["kind"];
	
	if (   estJson.find("forwardWalkingSteps") != estJson.end() )
	  {
	    int i=0;
	    
	    for ( auto &   fwStepJ : estJson["forwardWalkingSteps"]  )
	      {
		json_t jFW = estJson;
		int steps = fwStepJ.get<int>();
		
		jFW["forwardWalkingSteps"]=steps;
		jFW["targetLabel"]=jFW["label"];
		jFW["label"]=jFW["label"].get<std::string>() + "_fw" + std::to_string(steps);

		if (knownObservableTypes.at(id) == "scalar")
		  {
		    estimators.push_back(new realScalarForwardWalkingEstimator(jFW) );
		  }
		else if (knownObservableTypes.at(id) == "histogram")
		  {
		    estimators.push_back(new realHistogramForwardWalkingEstimator(jFW) );
		  }
		else
		  {
		    throw missingImplementation("Forward walking estimator for " + knownObservableTypes.at(id) + " not yet implemented");
		    
		  }
		i++;
	      }

	   
	    
	  }
	else
	  {
	    
	
	    if ( (id != "forceEnergy") and (id !="energy") )
	      {
		
		estimators.push_back( abstractFactory_t::create(id,estJson) );
		
		
	      }
	  }
      }
    
    return estimators;
  }
