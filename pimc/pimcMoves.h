#include "../src/traits.h"



class centerOfMassMove
{
    public:
        centerOfMassMove(real_t delta); // performs a move of the center of mass
        void move(const configurations_t & oldConf, configurations_t & newConf); // actually creates the new configurations from the old ones
        void acceptanceRatio(configurations_t & oldConf, configurations_t & newConf , action_t  & currentAction); // acceptance ratio of the current move
        
    private:

};

