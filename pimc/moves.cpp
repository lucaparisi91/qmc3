#include "moves.h"
#include "action.h"
#include "qmcExceptions.h"


namespace pimc
{      


    std::array<int,2> timeSliceGenerator::operator()(randomGenerator_t & randG, int nBeads , int maxBeadLength)
    {
        int t0=std::floor( uniformRealNumber(randG)*nBeads );
        int length = uniformRealNumber(randG)*maxBeadLength;

        int t1 = t0 + length;

        return {t0,t1};
    }


    void levyReconstructor::apply (configurations_t & configurationsNew, configurations_t & configurationsOld, std::array<int,2> timeRange,int iChain ,Real timeStep,randomGenerator_t & randG)
    {
        int l = timeRange[1] - timeRange[0];

        const auto & dataOld = configurationsOld.dataTensor() ;
        auto & dataNew = configurationsNew.dataTensor() ;
        
        const auto & chainsInfo = configurationsNew.getChainsInfo();

        int iChainNext = chainsInfo[iChain].getNextChain().index();

        auto timeRanges = splitPeriodicTimeSlice(timeRange,configurationsOld.nBeads());


        // copy initial and final configurations to the buffer
        if ( timeRange[1] >= configurationsOld.nBeads() )
        {
        
            for (int d=0;d<getDimensions();d++)
            {
                buffer(l,d)=dataOld(iChainNext,d,timeRanges[1][1]);
            }

        }
        else
        {
            for (int d=0;d<getDimensions();d++)
            {
                buffer(l,d)=dataOld(iChainNext,d,timeRange[1]);
            }

        }

        for (int d=0;d<getDimensions();d++)
            {
                buffer(0,d)=dataOld(iChain,d,timeRange[0]);
            }


        for (int t=0;t<l-1;t++)
        {
           
            Real eta = gauss(randG);

            for (int d=0;d<getDimensions();d++)
                {
                mean[d] = 
                (
                buffer(t ,d )*(l-t-1) 
                + buffer(l,d ) )
                /( l - t) 
                ;

                Real variance = (l-t-1) * 1. /(l-t) *timeStep;

                buffer( t + 1,d ) = mean[d] + eta *sqrt(variance) ;
                }
        }

        // copy the generated data back to the new configurations
        configurationsNew.copyDataFromBuffer( buffer, timeRanges[0], iChain   );
        configurationsNew.copyDataFromBuffer( buffer, timeRanges[1], iChain , timeRanges[0][1] - timeRanges[0][0] );

    }

levyMove::levyMove(levyReconstructor & levy_, int maxBeadLength_) : _levy(levy_) , uniformRealNumber(0,1),maxBeadLength(maxBeadLength_) , tmp(maxBeadLength_,getDimensions() ) {}



bool levyMove::attemptMove( configurations_t & confs, firstOrderAction & ST,randomGenerator_t & randG)
{

    int nChains = confs.nChains();
    int nBeads = confs.nBeads();
    auto & S = ST.getPotentialAction();

    int iChain=confsSampler.sampleChain(confs,randG);

    auto timeRange = tGen(randG,nBeads,maxBeadLength);

    auto timeRanges = splitPeriodicTimeSlice(timeRange,confs.nBeads());
    const auto & currentChain = confs.getChainsInfo()[iChain];

    int iChainNext = currentChain.getNextChain().index();

    if (not currentChain.contains(timeRange) ) // if crosses a head or tail
    {
        return false;
    }
    
    
    auto & data = confs.dataTensor();

    auto sOld = S.evaluate(confs, timeRanges[0], iChain);
    sOld += S.evaluate(confs, timeRanges[1], iChainNext);


    // copy to internal buffer beads to move
    copyToBuffer(confs,timeRanges[0],iChain);
    copyToBuffer(confs,timeRanges[1],iChainNext, timeRanges[0][1]  - timeRanges[0][0] );


    _levy.apply(confs,timeRange,iChain,S.getTimeStep(),randG);

    auto  sNew= S.evaluate(confs,timeRanges[0], iChain) ;
    sNew+=S.evaluate(confs,timeRanges[1], iChainNext) ;

    const auto actionDifference = sNew - sOld;

    bool accepted = sampler.acceptLog(-actionDifference,randG);

    if (! accepted)
    {
        // copy back old beads
        copyFromBuffer(confs,timeRanges[0],iChain);
        copyFromBuffer(confs,timeRanges[1],iChainNext, timeRanges[0][1]  - timeRanges[0][0] );
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
     auto & data = confs.dataTensor();
    // selects a worm at random
    auto & worms = confs.worms();
    int iWorm = std::floor(  uniformRealNumber(randG) * worms.size() );
    auto  oldWorm = worms[iWorm];

    auto & headChain = oldWorm.getHead();


    auto & Spot = S.getPotentialAction();

    // selects a time slice length at random. Refuse if reconstructed time slice crosses over the end bead

    int l = std::floor(uniformRealNumber(randG) * maxStepLength);

    std::array<int,2> timeSlice = {headChain.getHead(),headChain.getHead() + l};

    if (timeSlice[1] >= confs.nBeads() )
    {
        return false;
    } // no time wrap allowed

    
    // tower sampling a particle i with gaussian weights on the relative distances

   
    const auto & group = confs.getGroup(headChain.index() );

    std::array<Real, 3> distance;
    chainSampler.reset();

    for(int i=group.iStart;i<=group.iEnd;i++)
    {
        Real norm=0;
        for(int d=0;d<getDimensions();d++)
        {
            distance[d]=geo.difference(  data(headChain.index(),d,timeSlice[0]) - data(i,d,timeSlice[1]),d);
        }

        particleSelectionWeight=exp(freeParticleLogProbability(distance,S.getTimeStep(),group.mass));
        chainSampler.accumulateWeight(particleSelectionWeight);
    }

    int iPartner=particleSampler.sample(randG);

    const auto & chainsInfo = confs.getChainsInfo();

    const auto & partnerChain = chainsInfo[iPartner];

    if ( partnerChain.isOpen() )
    {
        return false;
    }

    Real deltaS=0;
    deltaS-=Spot.evaluate(confs,timeSlice,iPartner);

    // perform levy reconstruction
    for(int d=0;d<getDimensions();d++)
        {
            data(iPartner,d,timeSlice[1])=data(headChain.index(),d,timeSlice[1]);
        }

    _levy.apply(confs,  timeSlice,iPartner, Spot.getTimeStep() ,randG );


    deltaS-=Spot.evaluate(confs,timeSlice,iPartner,Spot.getTimeStep() );


    bool accept = metropolisSampler.acceptLog(-deltaS,randG);
    if (accept) 
    {
        // swap partner and head chain data below head time 
        confs.swapData({ headChain.getTail() + 1 ,headChain.getHead()-1},iPartner,headChain.getHead() );

        // relink swapped chains
        int iPrevPartner = partnerChain.getPrevChain().index();
        confs.join(headChain.getPrevChain().index(),iPartner);
        confs.join(iPrevPartner,headChain.index() );

    }
    else
    {
        
    }    

    return accept;
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


bool openMove::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{

    int iChain = confsSampler.sampleChain(confs,randG);
    Real var=2*D*timeStep;

    int iTime= std::floor( uniformRealNumber(randG) * confs.nBeads() );

    const auto & headChain = confs.getChainsInfo()[iChain];

    
    // creates a worm
    confs.open(iTime,iChain);
    
    const auto & currentWorm=confs.worms()[0];
    
    int iChainHead=currentWorm.getHead().index();
    int iChainTail=currentWorm.getTail().index();

    auto & data = confs.dataTensor();
    Real distanceSquared=0;
    const auto & geo = S.getGeometry();
    

    std::array<Real,3> difference;
    for (int d=0;d<getDimensions();d++)
    {
        data(iChainTail,d,iTime)=gauss(randG) + data(iChainHead,d,iTime);
        difference[d]=
            geo.difference( 
                data(iChainHead,d,iTime)-data(iChainTail,d,iTime),d
            );
    }


    action & potS = S.getPotentialAction();

    std::array<int,2> timeSlice={iTime,iTime};
    auto deltaS=0; // at first first order there is no change in potential action

    Real mass = 1.;


    auto propRatio = -deltaS - freeParticleLogProbability(difference,S.getTimeStep(),mass);

    bool accept = sampler.acceptLog(propRatio,randG);

    if ( accept)
    {
        confs.copyData({iTime+1,confs.nBeads()-1},iChainHead,iChainTail);
    }
    else
    {
        confs.close(0);
    }

    return accept;

};

bool closeMove::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{
    
    int iWorm = 0;
    Real var=2*D*timeStep;
    const auto & currentWorm=confs.worms()[0];

    int iChainHead=currentWorm.getHead().index();
    int iChainTail=currentWorm.getTail().index();


    int iTime=currentWorm.getHead().getHead();

    if ( iTime != currentWorm.getTail().getTail() )
    {
        throw missingImplementation("Head and tail have to be the same in the close move fpr now ");
    }


    auto & data = confs.dataTensor();
    const auto & geo = S.getGeometry();
    

    std::array<Real,3> difference;
    for (int d=0;d<getDimensions();d++)
    {
        data(iChainTail,d,iTime)=gauss(randG) + data(iChainHead,d,iTime);
        difference[d]=
            geo.difference( 
                data(iChainHead,d,iTime)-data(iChainTail,d,iTime),d
            );
    }


    action & potS = S.getPotentialAction();

    std::array<int,2> timeSlice={iTime,iTime};
    auto deltaS=0; // at first first order there is no change in potential action

    Real mass = 1.;


    auto propRatio = - deltaS +freeParticleLogProbability(difference,S.getTimeStep(),mass);

    bool accept = sampler.acceptLog(propRatio,randG);

    if ( accept)
    {
        confs.copyData({iTime+1,confs.nBeads()-1},iChainTail,iChainHead);
    }
    else
    {
        confs.close(0);
    }

    return accept;
    
};


advanceRecedeMove::advanceRecedeMove(int maxAdvanceLength_) :
maxAdvanceLength(maxAdvanceLength_),buffer(maxAdvanceLength_,getDimensions()) , _levy(maxAdvanceLength)
{

}

bool advanceRecedeMove::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG)
{
    // select a worm and generate a time slice to be reconstructed
    if (confs.worms().size() ==0 )
    {
        return false;
    }

    int iWorm= std::floor(distr(randG) * confs.worms().size() );
    auto & currentWorm=confs.worms()[iWorm];

    int iChainHead=currentWorm.getHead().index();
    int iChainTail=currentWorm.getTail().index();


    int l = std::floor(distr(randG)*maxAdvanceLength);

    std::array<int,2> timeRangeHead = {currentWorm.getHead().getHead(),currentWorm.getHead().getHead() + l};
    std::array<int,2> timeRangeTail = {currentWorm.getTail().getTail(),currentWorm.getTail().getTail() + l};
    
    

    auto timeRangesHead = splitPeriodicTimeSlice( timeRangeHead,confs.nBeads());
    auto timeRangesTail = splitPeriodicTimeSlice( timeRangeTail,confs.nBeads());


    Real deltaS=0;

    auto & sPot= S.getPotentialAction();

    // action of tail to be masked
    deltaS-=sPot.evaluate(confs,timeRangesTail[0],iChainTail);
    deltaS-=sPot.evaluate(confs,timeRangesTail[1],iChainTail);

    // copy the head bead
    auto & data = confs.dataTensor();
    int tBeadHead= (timeRangeHead[1] >= confs.nBeads() ) ? timeRangeHead[1] - confs.nBeads() : timeRangeHead[1];

    int iHeadBead= (timeRangeHead[1] >= confs.nBeads() ) ? iChainHead : iChainTail;

    for(int d=0;d<getDimensions();d++)
    {
        tmpMean[d]=data(iHeadBead,d,tBeadHead);
    }

    confSampler.sampleFreeParticlePosition(tmpPosition,tmpMean, S.getTimeStep() * l, randG  );

    for(int d=0;d<getDimensions();d++)
    {
        data(iChainHead,d,tBeadHead)=tmpPosition[d];
    }

    _levy.apply(confs,timeRangeHead,iChainHead,sPot.getTimeStep(),randG);


    // update worm configuration
    confs.moveTail(iChainHead,l);
    confs.moveTail(iChainTail,l);

    // evaluate potential due to reconstructed head beads
    deltaS+=sPot.evaluate(confs,timeRangesHead[0],iChainHead);
    deltaS+=sPot.evaluate(confs,timeRangesHead[1],iChainTail);

    bool accept = sampler.acceptLog(-deltaS,randG);

    if (not accept)
    {
        confs.moveHead(iWorm,-l);
        confs.moveTail(iWorm,-l);    
    }

    return accept;

}


}