#include "action.h"
namespace pimc
{

Real kineticAction::evaluate( pimcConfigurations_t & configurations , std::array<int,2> timeSlices , std::array<int,2> chainRange  )
{

    const auto & data = configurations.dataTensor();

    geo.updateSpringDifferences(distancesBuffer,data , timeSlices , chainRange );

    auto sum= reduceOnSpringDistances( [] (Real x,Real y,Real z){ return x*x + y*y + z*z;} ,distancesBuffer,timeSlices, chainRange,configurations.getMask());


    return sum/(4* D * tau);
    
}

Real kineticAction::evaluate( pimcConfigurations_t & configurations , std::array<int,2> timeSlices , int iChange )
{   
    Real sum=0;

    for (const auto & group : _particleGroups)
    {
        if ( ( iChange >= group.iStart )and iChange <= group.iEnd  )
        {
            sum+=group.mass*evaluate(configurations,timeSlices, {iChange,iChange});
        }
    }
    return sum;
}

Real kineticAction::evaluate( pimcConfigurations_t & configurations )
{
    Real sum=0;
    for (auto & group : _particleGroups)
    {
        sum+=evaluate(configurations  , {0 , nBeads-1} , {group.iStart,group.iEnd} ) ;
    }
    return sum;
};

}