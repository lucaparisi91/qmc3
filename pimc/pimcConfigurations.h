
#include <vector>
#include "../src/traits.h"

class pimcConfigurations
{
    public:
    pimcConfigurations(size_t NBeads, std::vector<int> nParticles ) ;
    size_t nBeads() const ;
    updateDistances(int iSet, int iParticle)
    
    private:
    std::vector< states_t  > positions;
    bool isOpenSector ;
    int iWorm,  headWorm ,  tailWorm;
    std::vector< tableDistances_t   > sameTimeDiffTable; // vector of distances at the same imaginary time slime
    std::vector<tableDistances_t >  sameParticleDiffTable; // distances between particles at different number but same chain


}