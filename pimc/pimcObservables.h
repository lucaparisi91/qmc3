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


class virialEnergyEstimator
{
    public:
    virialEnergyEstimator(int nMax, int MMax) : buffer(nMax,getDimensions(),MMax),rC(nMax,getDimensions(),MMax) {}
    Real operator()(configurations_t & configurations, firstOrderAction & S);
    
    private:
    Eigen::Tensor<Real,3> buffer;
    Eigen::Tensor<Real,3> rC;
};




}