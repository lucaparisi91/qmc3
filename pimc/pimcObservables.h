#include "tools.h"
#include "pimcConfigurations.h"
#include "action.h"


namespace pimc
{

class thermodynamicEnergyEstimator
{
    public:
    thermodynamicEnergyEstimator(){}
    Real operator()(configurations_t & configurations, firstOrderAction & S);
};

}