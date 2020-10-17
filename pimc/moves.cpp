#include "moves.h"
#include "action.h"

namespace pimc
{
    

    std::array<std::array<int,2>, 2> timeSliceGenerator::operator()(randomGenerator_t & randG, int nBeads , int maxBeadLength)
    {
        int t0=std::floor( uniformRealNumber(randG)*nBeads );
        int length = uniformRealNumber(randG)*maxBeadLength;

        int t1 = t0 + length;
        std::array<std::array<int,2>, 2> timeSlice;

        if ( t1 < nBeads )
        {
            timeSlice[0]={t0,t1};
            timeSlice[1]={0,-1};
        }
        else 
        {
            timeSlice[0]={t0, nBeads-1   };
            timeSlice[1]={0,t1%nBeads};
        }

        return timeSlice;

    }

    void levyReconstructor::apply (configurations_t & configurationsNew, configurations_t & configurationsOld, int iChain , std::array<int,2> timeRange)
    {
        int l = timeRange[1] - timeRange[0];
        const auto & dataOld = configurationsOld.dataTensor() ;
        auto & dataNew = configurationsNew.dataTensor() ;
        auto T=configurationsNew.nBeads();


        for (int d=0;d<getDimensions();d++)
        {
            dataNew(iChain,d,timeRange[0]%T)=dataOld(iChain,d,timeRange[0]%T);
            dataNew(iChain,d,timeRange[1]%T)=dataOld(iChain,d,timeRange[1]%T);
        }
        for (int t=0;t<l-1;t++)
        {
           
            Real eta = gauss(randG);

            for (int d=0;d<getDimensions();d++)
                {
                mean[d] = dataNew(iChain ,d, ( timeRange[0] + t )%T ) + (l  - t) *dataOld(iChain,d , timeRange[1]%T)/( l - t) ;
                Real variance = (l-t-1) /(l-t) *timeStep;

                dataNew(iChain,d, (timeRange[0] +t + 1)%T) = mean[d] + eta *sqrt(variance) ;
                }
        }

    }


levyMove::levyMove(levyReconstructor & levy_, int maxBeadLength_) : _levy(levy_) , uniformRealNumber(0,1),maxBeadLength(maxBeadLength_) {}


bool levyMove::attemptMove( configurations_t & confs, firstOrderAction & S)
{

    int nChains = confs.nChains();
    int nBeads = confs.nBeads();



    int iChain = std::floor( uniformRealNumber(randG) * nChains );


    std::array<int,2> timeRange;

    auto timeRanges = tGen(randG,nBeads,maxBeadLength);



    auto & data = confs.dataTensor();

    const auto sOld = S.evaluate(confs,timeRanges,iChain);

    // copy to internal buffer beads to move
    copyToBuffer(confs,timeRanges,iChain);


    const auto  sNew= S.evaluate(confs,timeRange, iChain) ;

    const auto actionDifference = sNew - sOld;

    bool accepted = sampler.acceptLog(actionDifference,randG);

    if (! accepted)
    {
        // copy old beads from initial data
        copyToBuffer(confs,timeRanges, iChain);

    }

    return accepted;

}

void levyMove::copyToBuffer(configurations_t & confs, std::array<int,2> timeRange,int iChain)
{
    auto & data = confs.dataTensor();
    tmp.resize(confs.nBeads(),getDimensions()) ;

    
    for (int t=timeRange[0] ;t<=timeRange[1];t++ )
        {
            for (int d=0;d<getDimensions();d++)
            {
            tmp(t,d)=data(iChain,d,t);
            }
        }

};

void levyMove::copyFromBuffer(configurations_t & confs, std::array<int,2> timeRange,int iChain)
{
    auto & data = confs.dataTensor();


    for (int t=timeRange[0] ;t<=timeRange[1];t++ )
        {
            for (int d=0;d<getDimensions();d++)
            {
            data(iChain,d,t)=tmp(t,d);
            }
        }
};


}