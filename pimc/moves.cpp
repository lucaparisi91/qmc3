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



    void levyReconstructor::apply (configurations_t & configurationsNew, configurations_t & configurationsOld, std::array<int,2> timeRange,int iChain ,Real timeStep,randomGenerator_t & randG)
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

bool levyMove::isValidSlice( configurations_t & confs , const std::array<int,2> timeRange,int iChain) const
{
    int iWorm=confs.findWorm(iChain);
    if (iWorm >0)
    {
        auto & worm= confs.worms()[iWorm];
        if (worm.iChainHead == iChain )
        {
            if(timeRange[1]>worm.iHead)
                {
                    return false;
                }
            else
            {
                return true;
            }
            
        }   
        else if (worm.iChainTail == iChain )
        {
            if ( 
                (timeRange[0]<worm.iTail) or 
                (timeRange[1]>(confs.nBeads()-1) )
                )
                {
                    return false;
                }
                else
                {
                    return true;
                }
        }

        else
        {
            throw invalidState("Chain does not seem to be a worm.");
        }

    }

    return true;    
}


bool levyMove::attemptMove( configurations_t & confs, firstOrderAction & ST,randomGenerator_t & randG)
{

    int nChains = confs.nChains();
    int nBeads = confs.nBeads();
    auto & S = ST.getPotentialAction();
    
    int iChain=confsSampler.sampleChain(confs,randG);


    auto timeRange = tGen(randG,nBeads,maxBeadLength);

    if (not isValidSlice(confs,timeRange,iChain) )
    {
        return false;
    }

    auto timeRanges = splitPeriodicTimeSlice(timeRange,nBeads);
    
    auto & data = confs.dataTensor();

    const auto sOld = S.evaluate(confs,timeRanges,iChain);

    // copy to internal buffer beads to move
    copyToBuffer(confs,timeRanges,iChain);

    _levy.apply(confs,timeRange,iChain,S.getTimeStep(),randG);

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
     auto & data = confs.dataTensor();
    // selects a worm at random
    auto & worms = confs.worms();
    int iWorm = std::floor(  uniformRealNumber(randG) * worms.size() );
    auto  oldWorm = worms[iWorm];

    auto & Spot = S.getPotentialAction();

    // selects a time slice length at random. Refuse if reconstructed time slice crosses over the end bead


    int l = std::floor(uniformRealNumber(randG) * maxStepLength);

    std::array<int,2> timeSlice = {oldWorm.iHead,oldWorm.iHead + l};

    if (timeSlice[1] >= confs.nBeads() )
    {
        return false;
    } // no time wrap allowed

    
    // tower sampling a particle i with gaussian weights on the relative distances

   
    const auto & group = confs.getGroup(oldWorm.iChainHead);

    std::array<Real, 3> distance;
    chainSampler.reset();

    for(int i=group.iStart;i<=group.iEnd;i++)
    {
        Real norm=0;
        for(int d=0;d<getDimensions();d++)
        {
            distance[d]=geo.difference(  data(oldWorm.iChainHead,d,timeSlice[0]) - data(i,d,timeSlice[1]),d);
        }

        particleSelectionWeight=exp(freeParticleLogProbability(distance,S.getTimeStep(),group.mass));
        chainSampler.accumulateWeight(particleSelectionWeight);
    }

    int iPartner=particleSampler.sample(randG);
    if (confs.isWorm(iPartner) )
    {
        return false;
    }

    Real deltaS=0;
    deltaS-=Spot.evaluate(confs,timeSlice,iPartner);

    // perform levy reconstruction
    for(int d=0;d<getDimensions();d++)
        {
            data(oldWorm.iChainHead,d,timeSlice[1])=data(iPartner,d,timeSlice[1]);
        }

    _levy.apply(confs,  timeSlice,oldWorm.iChainHead , Spot.getTimeStep() ,randG );

    //updates the worm
    worm newWorm(oldWorm);
    newWorm.iChainHead=iPartner;
    confs.updateWorm(iWorm,newWorm);

    deltaS-=Spot.evaluate(confs,timeSlice,oldWorm.iChainHead,Spot.getTimeStep() );


    bool accept = metropolisSampler.acceptLog(-deltaS,randG);

    if (accept)
    {
        confs.copyData({timeSlice[1],confs.nBeads()-1},iPartner,oldWorm.iChainHead);
    }
    else
    {
        confs.updateWorm(iWorm,oldWorm);
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

    // creates the additional bead
    
    int iWorm=confs.open(iTime,iChain,false);
    auto worm = confs.worms()[iWorm];
    auto & data = confs.dataTensor();
    Real distanceSquared=0;

    for (int d=0;d<getDimensions();d++)
    {
        data(worm.iChainHead,d,iTime)=gauss(randG) + data(iChain,d,iTime);
        distanceSquared+=std::pow( data(worm.iChainHead,d,iTime)-data(worm.iChainTail,d,iTime),2);
    }



    action & potS = S.getPotentialAction();

    std::array<int,2> timeSlice={iTime,iTime};
    auto deltaS=potS.evaluate(confs,timeSlice,worm.iChainHead);

    auto propRatio = -deltaS + 0.5*log(2*M_PI*var) + 0.5*distanceSquared/(var) + log(C); // proposed ratio for the new proposal

    bool accept = sampler.acceptLog(propRatio,randG);

    if ( accept)
    {
        // copy  the lower half of the old chain
        for(int t=0; t< iTime ; t++)
        {
            for(int d=0;d<getDimensions();d++)
            {
                  data(worm.iChainHead,d,t)=data(worm.iChainTail,d,t);
            }
          
        }
    }
    else
    {
        confs.close(iWorm,false);
    }

    return accept;

};


bool closeMove::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{
    
    int iWorm = 0;
    Real var=2*D*timeStep;
    auto worm = confs.worms()[iWorm];

    auto & data = confs.dataTensor();
    Real distanceSquared=0;

    for (int d=0;d<getDimensions();d++)
    {

        distanceSquared+=std::pow( data(worm.iChainHead,d,worm.iHead)-data(worm.iChainTail,d,worm.iTail),2);
    }

    action & potS = S.getPotentialAction();

    std::array<int,2> timeSlice={worm.iHead,worm.iHead};
    auto deltaS=-potS.evaluate(confs,timeSlice,worm.iChainHead);

    auto propRatio = -deltaS - 0.5*log(2*M_PI*var) - 0.5*distanceSquared/(var) - log(C) ; // proposed ratio for the new proposal

    bool accept = sampler.acceptLog(propRatio,randG);

    if ( accept)
    {
        confs.close(iWorm,true);
    }

    return accept;
};


advanceRecedeMove::advanceRecedeMove(int maxAdvanceLength_) :
maxAdvanceLength(maxAdvanceLength_),buffer(maxAdvanceLength_,getDimensions())
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
    auto oldWorm=confs.worms()[iWorm];

    int l = std::floor(distr(randG)*maxAdvanceLength);

    std::array<int,2> timeRange = {oldWorm.iHead,oldWorm.iHead + l};
    auto timeRanges = splitPeriodicTimeSlice(timeRange,confs.nBeads());

    // copy head to be overwritten to buffer. Only if time slice wraps around the time boundary    
    confs.copyDataToBuffer(buffer,timeRanges[1],oldWorm.iChainHead);

    Real deltaS=0;

    
     auto & sPot= S.getPotentialAction();

    // compute the potential action of tail and head to be overwritten    
    std::array<int,2> tailSlice={confs.nBeads()- (timeRanges[1][1]-timeRanges[1][0] + 1),confs.nBeads()-1};
    deltaS-=sPot.evaluate(confs,tailSlice,oldWorm.iChainTail);
    deltaS-=sPot.evaluate(confs,timeRanges[1],oldWorm.iChainHead);
    
    //configurations_t::copyData(confs,tailSlice,{oldWorm.iChainTail,oldWorm.iChainTail},confs,timeRanges[1][0] ,oldWorm.iChainHead);

    
    auto & data = confs.dataTensor();

    int tHead=timeRange[1]%confs.nBeads();

    for(int d=0;d<getDimensions();d++)
    {
        tmpMean[d]=data(oldWorm.iChainHead,d,tHead);
    }

    confSampler.sampleFreeParticlePosition(tmpPosition,tmpMean, S.getTimeStep() * l, randG  );

    for(int d=0;d<getDimensions();d++)
    {
        data(oldWorm.iChainHead,d,tHead)=tmpPosition[d];
    }

    _levy.apply(confs,timeRange,oldWorm.iChainHead,sPot.getTimeStep(),randG);

    //copy the overwirtten beads at the beginning head to the beginning of the tail
    configurations_t::copyData(confs,timeRanges[1],{oldWorm.iChainHead,oldWorm.iChainHead},confs,0,oldWorm.iChainTail);

    // update worm configuration
    worm newWorm(oldWorm);

    if (timeRange[1]< confs.nBeads()) // does not wrap time boundary
    {
        newWorm.iHead=timeRange[1];
        newWorm.iTail=timeRange[1];
    }
    else
    {
        newWorm.iHead=timeRanges[0][1];
        newWorm.iTail=timeRanges[0][1];
        std::swap(newWorm.iChainHead,newWorm.iChainTail);

    }
    
    confs.updateWorm(iWorm,newWorm);

    // evaluate potential due to reconstructed head beads. If wraps time boundary some will be in the tail
    deltaS+=sPot.evaluate(confs,timeRanges[0],oldWorm.iChainHead);
    deltaS+=sPot.evaluate(confs,timeRanges[1],oldWorm.iChainTail);

    bool accept = sampler.acceptLog(-deltaS,randG);

    if (not accept)
    {
        confs.updateWorm(iWorm,oldWorm);
        confs.copyDataFromBuffer(buffer,timeRanges[1],oldWorm.iChainHead);
    }


    return accept;






}


}