#include "toolsPimc.h"

namespace pimc
{

std::array<std::array<int,2>, 2> splitPeriodicTimeSlice(const std::array<int,2> & timeSlice, int nBeads)
    {
        auto [t0,t1] = timeSlice;
        std::array<std::array<int,2>, 2> timeSlices;

        if ( t1 < nBeads )
        {
            timeSlices[0]={t0,t1};
            timeSlices[1]={0,-1};
        }
        else 
        {
            timeSlices[0]={t0, nBeads-1   };
            timeSlices[1]={0,t1%nBeads};
        }

        return timeSlices;

    }



Real freeParticleLogProbability(std::array<Real,3> & delta,Real tau,Real mass)
    {
        const Real D = 0.5;
        Real var = 2 * D * tau / mass;

        Real p=0;
        
        for(int d=0;d<getDimensions();d++)
        {
            p+= delta[d]*delta[d];
        }
        p*=-0.5 /var;
        p+= -0.5*log(2*M_PI*var);
        
        return p;
    }


}