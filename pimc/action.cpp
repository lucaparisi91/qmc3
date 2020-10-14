#include "action.h"
namespace pimc
{

Real kineticAction::evaluate( pimcConfigurations_t & configurations , std::array<int,2> timeSlices , std::array<int,2> chainRange  )
{

    const auto & data = configurations.dataTensor();

    geo.updateSpringDifferences(distancesBuffer,data , timeSlices , chainRange );

    auto sum= reduceOnSpringDistances( [] (Real x,Real y,Real z){ return x*x + y*y + z*z;} ,distancesBuffer,timeSlices, chainRange);

    return 0;

}

Real kineticAction::evaluate( pimcConfigurations_t & configurations )
{
    auto sum=evaluate(configurations  , {0 , nBeads} , {0,nChains-1} ) ;

    return sum;
}



}