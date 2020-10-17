#include "pimcConfigurations.h"
namespace pimc
{

void pimcConfigurations::deleteBeads(   std::array<int,2> timeRange, int iChain )
{
    for (int t=timeRange[0] ; t<= timeRange[1] ;t++   )
    {
        _mask(t,iChain)=0;
    }
}

void pimcConfigurations::createBeads(   std::array<int,2> timeRange, int iChain )
{
    for (int t=timeRange[0] ; t<= timeRange[1] ;t++   )
    {
        _mask(t,iChain)=1;
    }
}




};

