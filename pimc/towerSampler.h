#include "tools.h"

namespace pimc
{
class towerSampler
{
    public:
    towerSampler() : uniformRealNumber(0,1),totWeight(0) {} 
    
    towerSampler(int max_num_weights) : accWeights(max_num_weights) {}


    virtual void reset() {totWeight=0;int iCurrentWeight;}

    int sample(std::vector<Real> & accWeights, Real totWeight, randomGenerator_t & randG)
    {
        Real randomWeight=uniformRealNumber(randG) * totWeight;
        int i=0;
        for(i=0; accWeights[i]<randomWeight;i++) {}

        return i;
    }

    private:
    std::uniform_real_distribution<float> uniformRealNumber;
    std::vector<Real> accWeights;
    Real totWeight;

};

}