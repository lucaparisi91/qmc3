#include "action.h"

namespace pimc
{

Real kineticAction::evaluate( pimcConfigurations_t & configurations , std::array<int,2> timeSlices , std::array<int,2> chainRange  )
{
    const auto & data = configurations.dataTensor();
    auto & geo=getGeometry();

    geo.updateSpringDifferences(distancesBuffer,data , timeSlices , chainRange );
    #if DIMENSIONS == 3
    auto sum= reduceOnSpringDistances( [] (Real x,Real y,Real z){ return x*x + y*y + z*z;} ,distancesBuffer,timeSlices, chainRange,configurations.getMask());
    #endif
    #if DIMENSIONS == 2
    auto sum= reduceOnSpringDistances( [] (Real x,Real y){ return x*x + y*y ;} ,distancesBuffer,timeSlices, chainRange,configurations.getMask());
    #endif
    #if DIMENSIONS == 1
    auto sum= reduceOnSpringDistances( [] (Real x){ return x*x;} ,distancesBuffer,timeSlices, chainRange,configurations.getMask());
    #endif

    return sum/(4* D * getTimeStep());
}

Real kineticAction::evaluate( pimcConfigurations_t & configurations , std::array<int,2> timeSlices , int iChange )
{   
    Real sum=0;

    const auto & particleGroups = configurations.getGroups();

    for (const auto & group : particleGroups)
    {
        if ( ( iChange >= group.iStart )and iChange <= group.iEnd  )
        {
            sum+=group.mass*evaluate(configurations,timeSlices, {iChange,iChange});
        }
    }
    return sum;
}


Real kineticAction::evaluate( pimcConfigurations_t & configurations , std::array<int,2> timeSlices , int iChain1 , int iChain2 )
{   
    return evaluate(configurations,timeSlices,iChain1) + evaluate(configurations,timeSlices,iChain2) ; 
}


Real kineticAction::evaluate( pimcConfigurations_t & configurations )
{
    Real sum=0;
    
    const auto & particleGroups = configurations.getGroups();
    
    for (auto & group : particleGroups)
    {
        sum+=evaluate(configurations  , {0 , nBeads-1} , {group.iStart,group.iEnd} ) ;
    }

    return sum;
};


}