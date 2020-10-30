#include "moves.h"
#include "action.h"


namespace pimc
{      
    std::array<int,2> timeSliceGenerator::operator()(randomGenerator_t & randG, int nBeads , int maxBeadLength)
    {
        int t0=std::floor( uniformRealNumber(randG)*nBeads );
        int length = uniformRealNumber(randG)*maxBeadLength;

        int t1 = t0 + length;

        return {t0,t1};
    }

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


    void levyReconstructor::apply (configurations_t & configurationsNew, configurations_t & configurationsOld, int iChain , std::array<int,2> timeRange,randomGenerator_t & randG)
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
                mean[d] = 
                (
                dataNew(iChain ,d, ( timeRange[0] + t )%T )*(l-t-1) 
                + dataOld(iChain,d , timeRange[1]%T) )
                /( l - t) ;
                Real variance = (l-t-1) * 1. /(l-t) *timeStep;

                dataNew(iChain,d, (timeRange[0] +t + 1)%T) = mean[d] + eta *sqrt(variance) ;
                }
        }

    }


levyMove::levyMove(levyReconstructor & levy_, int maxBeadLength_) : _levy(levy_) , uniformRealNumber(0,1),maxBeadLength(maxBeadLength_) , tmp(maxBeadLength_,getDimensions() ) {}

bool levyMove::attemptMove( configurations_t & confs, firstOrderAction & ST,randomGenerator_t & randG)
{

    int nChains = confs.nChains();
    int nBeads = confs.nBeads();
    auto & S = ST.getPotentialAction();

    const auto & groups = confs.getGroups();
    int iGroup=std::floor(uniformRealNumber(randG)*groups.size());
    const auto & group = groups[iGroup];

    int iChain=std::floor(uniformRealNumber(randG)*group.size() ) +group.iStart;


    auto timeRange = tGen(randG,nBeads,maxBeadLength);

    auto timeRanges = splitPeriodicTimeSlice(timeRange,nBeads);
    
    auto & data = confs.dataTensor();

    const auto sOld = S.evaluate(confs,timeRanges,iChain);

    // copy to internal buffer beads to move
    copyToBuffer(confs,timeRanges,iChain);

    _levy.apply(confs,confs,iChain,timeRange,randG);

    const auto  sNew= S.evaluate(confs,timeRanges, iChain) ;

    const auto actionDifference = sNew - sOld;

    bool accepted = sampler.acceptLog(-actionDifference,randG);


    if (! accepted)
    {
        // copy old beads from initial data
        copyFromBuffer(confs,timeRanges, iChain);

    }

    return accepted;
}

void levyMove::copyToBuffer(configurations_t & confs, std::array<int,2> timeRange,int iChain, int offset)
{
    auto & data = confs.dataTensor();

    for (int t=timeRange[0] , tt=offset ;t<=timeRange[1];t++ & tt++ )
        {
            for (int d=0;d<getDimensions();d++)
            {
            tmp(tt,d)=data(iChain,d,t);
            }
        }
};

void levyMove::copyFromBuffer(configurations_t & confs, std::array<int,2> timeRange,int iChain,int offset)
{
    auto & data = confs.dataTensor();
    for (int t=timeRange[0],tt=offset ;t<=timeRange[1];t++ & tt++ )
        {
            for (int d=0;d<getDimensions();d++)
            {
                data(iChain,d,t)=tmp(tt,d);
            }
        }
};

bool swapMove::attemptMove(configurations_t & confs, firstOrderAction & S,randomGenerator_t & randG)
{
   /*  auto & data = confs.dataTensor();
    // selects a worm at random
    auto & worms = confs.worms();
    int iWorm = std::floor(  uniformRealNumber(randG) * worms.size() );
    auto  worm = worms[iWorm];

    auto & Spot = S.getPotentialAction();
    

    // selects a time slice length at random. Refuse if reconstructed time slice crosses over the end bead

    int timeSliceLength = std::floor(  uniformRealNumber(randG) * maxStepLength );
    int tEnd = worm.iHead + timeSliceLength; 
    if ( tEnd > confs.nBeads() ) return false;

    // tower sampling a particle i with gaussian weights on the relative distances


    int iWChain = worm.iChainHead;
    const auto & group = confs.getGroup(iWChain);


    int tStart=worm.iHead;
    particleSelectionWeight=0;

    Real timeStep = _levy.getTimeStep();
    for(int i=group.iStart;i<=group.iEnd;i++)
    {
        Real norm=0;
        for(int d=0;d<getDimensions();d++)
        {
            Real deltad=geo.difference(  data(iWChain,d,tStart) - data(i,d,tEnd),d);
            norm+=deltad*deltad;
        }

        particleSelectionWeight+=exp(-norm/(4*D*timeStep));
        particleSelectionAccWeights[i]=particleSelectionWeight;
    }

    int iPartner=particleSampler.sample(particleSelectionAccWeights, particleSelectionWeight, randG);

    if (iPartner == iWChain)
    {
        return false;
    }

    auto unChangedAction=Spot.evaluate(confs,{tStart,tEnd},iWChain,iWorm);


    // save particle partitions along the worm that will be ooverwritten

    for (int t=worm.iTail , tt=0;t<=tEnd ;t++ & tt++)
    {
        for(int d=0;d<getDimensions();d++)
        {
                tmpChainWorm(tt,d)=data(iWChain,d,t);
                tmpChainParticle(tt,d)=data(iPartner,d,t);
        }

    }

    // do levy reconstruction
    for(int d=0;d<getDimensions();d++)
    {
        data(iWChain,d,tEnd)=data(iPartner,d,tEnd);
    }

    for (int t=worm.iTail , tt=0;t<=tEnd ;t++ & tt++)
    {
        for(int d=0;d<getDimensions();d++)
        {
            data(iPartner,d,t)=tmpChainWorm(tt,d);
        }

    }

    _levy.apply(confs,confs, iWChain , {tStart,tEnd}, randG );

    confs.close(iWorm);
    int iNewWorm=confs.open(,iPartner);

     // compute the potential action on the beads
    auto changedAction=Spot.evaluate(confs,{tStart,tEnd},iWChain,iWorm);


    bool accept = metropolisSampler.acceptLog(- ( changedAction - unChangedAction),randG);


    if (accept)
    {
        // swap worm and partner particles above the thresold
        for(int t=worm.iTail+1;t<confs.nBeads();t++)
        {
            for (int d=0;d<getDimensions();d++)
            {
                std::swap(data(iWChain,d,t) , data(iPartner,d,t)    );
            }
        }
    }
    else
    {
        // revert to old old configuration
        confs.close(iNewWorm);
        confs.open({worm.iHead,worm.iTail},iWChain);

    }

    return accept;     */
}


void tableMoves::push_back(move * move_,Real weight)
{
    totalWeight+=weight;
    accumulatedWeights.push_back(totalWeight);
    _moves.push_back(move_);
    
};

move & tableMoves::sample(randomGenerator_t & randG)
{
    int iMove = sampler.sample(accumulatedWeights,totalWeight,randG);
    return *(_moves[iMove]);

};


}