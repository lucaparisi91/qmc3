
#include "toolsPimc.h"
#include "action.h"
#include "moves.h"

namespace pimc{

    class pimcDriver
    {
        public:
        pimcDriver(const json_t & j);
        void run();   


    private:
    int nBeads;
    firstOrderAction S;
    int stepsPerBlock;
    int nBlocks;
    int correlationSteps;
    std::vector<int> nParticles;
    int seed;
    Real timeStep;
    geometryPBC_PIMC geo;
    json_t j;
    
    tableMoves tab;

 };
}