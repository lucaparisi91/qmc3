#include "..src/traits/traits.h"


class action
{
    public:
    action(real_t tau);
    real_t operator()(configurations_t & conf);
    differenceKineticAction(int iSet, int iParticle , int iStartTime, int iEndTime, const configurations_t & oldConf,const configurations_t & newConf); // differences in the kinetic action when changing a strip of a worm
    differenceTwoBodyPotentialAction(int iSet, int iParticle , int iStartTime, int iEndTime, const configurations_t & oldConf,const configurations_t & newConf);

    private:
    potential_t _twoBodyPotential;
    real_t tau;
};

