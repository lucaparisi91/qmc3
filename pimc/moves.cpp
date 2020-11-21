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


    void levyReconstructor::apply (configurations_t & configurations, std::array<int,2> timeRange,int iChain ,Real timeStep,randomGenerator_t & randG)
    {
        /* Does not take into account time periodic boundary conditions
        Reconstruct between timeRange[1] -1 and timeRange[0] + 1
        timeRange[0] and timeRange[1] are used as boundary conditions
        */
        int l = timeRange[1] - timeRange[0];

        const auto & dataOld = configurations.dataTensor() ;
        auto & dataNew = configurations.dataTensor() ;
        
        auto timeRanges = splitPeriodicTimeSlice(timeRange,configurations.nBeads() );

        const auto & currentChain = configurations.getChain(iChain);
        int iChainNext=currentChain.next;



        for (int d=0;d<getDimensions();d++)
            {
                dataNew(iChain,d,timeRange[0])=dataOld(iChain,d,timeRange[0]);
            }


        // performs the actual copy
        for (int t=0;t<l-1;t++)
        {
           
            Real eta = gauss(randG);

            for (int d=0;d<getDimensions();d++)
                {
                mean[d] = 
                (
                dataNew(iChain,d,t+timeRange[0])*(l-t-1) 
                + dataNew(iChain,d , l + timeRange[0]  ) )
                /( l - t) 
                ;

                Real variance = (l-t-1) * 1. /(l-t) *timeStep;

                dataNew(iChain,d, timeRange[0] + t + 1 ) = mean[d] + eta *sqrt(variance) ;
                }
        }
        

        //configurationsNew.copyData(  {configurationsNew.nBeads() , timeRange[1]-1} , iChain , 0  , iChainNext);
        //configurationsNew.fillHead(iChain);

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
    const auto & currentChain = confs.getChain(iChain);

    int iChainNext = currentChain.next;

    if (
        (timeRanges[0][1] > currentChain.head )
         or (
            ( currentChain.next != -1 ) and
            ( confs.getChain(iChainNext).head <= timeRanges[1][1]   )
        )
    )
    {
        return false; // do not accept if the time range is not included inside a chain (including time PBC)


    }


    
    auto & data = confs.dataTensor();

    auto sOld = S.evaluate(confs, timeRanges[0], iChain);
    sOld += S.evaluate(confs, timeRanges[1], iChainNext);


    // copy to internal buffer beads to move
    copyToBuffer(confs,timeRanges[0],iChain);
    copyToBuffer(confs,timeRanges[1],iChainNext, timeRanges[0][1]  - timeRanges[0][0] + 1);


    _levy.apply(confs,timeRange,iChain,S.getTimeStep(),randG);
    confs.copyData( { timeRanges[0][1]+1 , timeRange[1] } , iChain, iChainNext ); // time periodic boundary conditions

    auto  sNew= S.evaluate(confs,timeRanges[0], iChain) ;
    sNew+=S.evaluate(confs,timeRanges[1], iChainNext) ;

    const auto actionDifference = sNew - sOld;

    bool accepted = sampler.acceptLog(-actionDifference,randG);

    if (! accepted)
    {
        // copy back old beads
        copyFromBuffer(confs,timeRanges[0],iChain);
        copyFromBuffer(confs,timeRanges[1],iChainNext, timeRanges[0][1]  - timeRanges[0][0] + 1 );
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

    int iChainHead = std::floor(  uniformRealNumber(randG) * confs.tails().size() );

    const auto  headChain=confs.getChain(iChainHead);



    auto & Spot = S.getPotentialAction();

    // selects a time slice length at random. Refuse if reconstructed time slice crosses over the end bead

    int l = std::floor(uniformRealNumber(randG) * maxStepLength);

    std::array<int,2> timeSlice = {headChain.head,headChain.head + l};

    if (timeSlice[1] >= confs.nBeads() )
    {
        return false;
    } // no time wrap allowed

    
    // tower sampling a particle i with gaussian weights on the relative distances
   
    const auto & group = confs.getGroup( iChainHead );

    std::array<Real, 3> distance;
    chainSampler.reset();

    for(int i=group.iStart;i<=group.iEnd;i++)
    {
        Real norm=0;
        for(int d=0;d<getDimensions();d++)
        {
            distance[d]=geo.difference(  data(iChainHead,d,timeSlice[0]) - data(i,d,timeSlice[1]),d);
        }

        particleSelectionWeight=exp(freeParticleLogProbability(distance,S.getTimeStep(),group.mass));
        chainSampler.accumulateWeight(particleSelectionWeight);
    }

    int iPartner=particleSampler.sample(randG);

    
    const auto  partnerChain = confs.getChain(iPartner);

    if ( partnerChain.isOpen() )
    {
        return false;
    }

    Real deltaS=0;
    deltaS-=Spot.evaluate(confs,timeSlice,iPartner);

    // perform levy reconstruction
    for(int d=0;d<getDimensions();d++)
        {
            data(iPartner,d,timeSlice[1])=data(iChainHead,d,timeSlice[1]);
        }

    _levy.apply(confs,  timeSlice,iPartner, Spot.getTimeStep() ,randG );


    deltaS-=Spot.evaluate(confs,timeSlice,iPartner,Spot.getTimeStep() );


    bool accept = metropolisSampler.acceptLog(-deltaS,randG);
    if (accept) 
    {
        // swap partner and head chain data below head time 
        confs.swapData({ headChain.tail + 1 ,headChain.head-1},iPartner,headChain.head );

        // relink swapped chains
        confs.swapTails(iChainHead,iPartner);

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

openMove::openMove(Real C_ , int maxReconstructedLength_) : C(C_), _levy(maxReconstructedLength_) ,  _maxReconstructedLength(maxReconstructedLength_) ,buffer(maxReconstructedLength_,getDimensions()){}


bool openMove::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{

    int iChain = confsSampler.sampleChain(confs,randG);
    Real var=2*D*timeStep;

    int iTime= std::floor( uniformRealNumber(randG) * confs.nBeads() );
    //int l= std::floor( uniformRealNumber(randG) * (confs.nBeads() -1) ) +1; // distance from itime where the head is formed

    int l=1;
    std::array<int,2> timeRange{iTime,iTime+l};


    

    //auto iTime=confs.nBeads();

    const auto headChain = confs.getChain(iChain);


    auto & data = confs.dataTensor();
    Real distanceSquared=0;
    const auto & geo = S.getGeometry();


    std::array<Real,3> difference;
    std::array<Real,3 > oldBead;

    // save the candidate head bead at time iTime + 1
    for (int d=0;d<getDimensions();d++)
    {
        oldBead[d]=data(iChain,d,iTime+1);
    }

    // generate the position of the new bead
    for (int d=0;d<getDimensions();d++)
    {
        
        data(iChain,d,iTime+1)=gauss(randG) + data(iChain,d,iTime+1);
        difference[d]=
            geo.difference( 
                data(iChain,d,iTime)-data(iChain,d,iTime+1),d
            );
    }

    auto deltaS=0;

    Real mass = confs.getGroupByChain(iChain).mass;

    auto propRatio = -deltaS - freeParticleLogProbability(difference,S.getTimeStep(),mass);

    bool accept = sampler.acceptLog(propRatio,randG);

    if ( accept)
    {

        if (iTime==confs.nBeads() - 1) // next chain will be the new tail
        {
            int tailChain=headChain.next;
            confs.setHead(iChain,iTime+1);
        }
        else // creates a new chain to contain the new tail
        {
            int tailChain=confs.pushChain(confs.getGroupByChain(iChain));
            confs.setTail(tailChain,iTime);
            confs.setHead(tailChain,confs.nBeads() );
            confs.setHead(iChain,iTime+1);
            confs.join(headChain.prev,tailChain);
            confs.copyData({iTime+1,confs.nBeads()-1}, iChain, tailChain );

            for (int d=0;d<getDimensions();d++) // copy back the original bead in the new tail
            {
                data(tailChain,d,iTime+1)=oldBead[d];
            }
            confs.fillHead(tailChain);
        }
    }
    else
    {
        for (int d=0;d<getDimensions();d++) // copy back the original bead
        {
            data(iChain,d,iTime+1)=oldBead[d]; 
        }

    }

    return accept;
};


closeMove::closeMove(Real C_ , int maxReconstructionLength) : C(C_),_levy(maxReconstructionLength),_maxLength(maxReconstructionLength),buffer(maxReconstructionLength,getDimensions()){}


bool closeMove::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{
    
    int iWorm = 0;
    Real var=2*D*timeStep;

    int iChainHead=std::floor( uniformRealNumber(randG)*confs.heads().size());
    int iChainTail=std::floor( uniformRealNumber(randG)*confs.tails().size());
    
    auto headChain=confs.getChain(iChainTail);
    auto tailChain=confs.getChain(iChainHead);

    auto & data = confs.dataTensor();
    const auto & geo = S.getGeometry();
   
    int iTime=headChain.head - 1;    

    int l = tailChain.tail - headChain.head + 1;

    l= l > confs.nBeads()/2. ? l=l-confs.nBeads() : 
    (l< -confs.nBeads()/2. ) ? l + confs.nBeads() : l;  // time pbc
    action & potS = S.getPotentialAction();
    Real deltaS=0;

    if (l< 0 )
    {
        throw missingImplementation("Head should be lower then the tail ");
    }
    else
    {
        std::array<int,2> timeRange={iTime,iTime + l};
        auto timeRanges= splitPeriodicTimeSlice(timeRange,confs.nBeads());

        // save configurations to buffer
        confs.copyDataToBuffer(buffer,timeRanges[0],iChainHead);
        confs.copyDataToBuffer(buffer,timeRanges[1],timeRanges[0][1] - timeRanges[0][0]+1);


        // apply the reconstruction
        for(int d=0;d<getDimensions();d++)
        {
            data(iChainHead,d,timeRange[1] )=data(iChainTail,d,tailChain.tail+1);
        }
        _levy.apply(confs, timeRange,iChainHead,S.getTimeStep(),randG);
        confs.copyData( {confs.nBeads() , timeRange[0]-1 } , iChainHead , 0,iChainTail  );

        // evaluate the difference in potential action
        deltaS-=potS.evaluate(confs,timeRanges[0],iChainHead);
        std::array<int,2> rangeWrapped=timeRanges[1];
        rangeWrapped[1]=-1;
        deltaS-=potS.evaluate(confs,rangeWrapped,iChainTail);

        // evaluates the free particle propagator
        std::array<Real,3> difference;
        for (int d=0;d<getDimensions();d++)
        {
            difference[d]=
                geo.difference( 
                    data(iChainHead,d,iTime)-data(iChainTail,d,iTime),d
                );
        }

        Real mass = 1.;
        auto propRatio = - deltaS +freeParticleLogProbability(difference,S.getTimeStep(),mass);
        bool accept = sampler.acceptLog(propRatio,randG);

        if ( accept )
        {
            // change positions of heads and tails
            if ( timeRange[1] >= confs.nBeads() )
            { // join the two chains
                confs.setHead(iChainHead,confs.nBeads() );
                confs.setTail(iChainTail,-1 );
                confs.join(iChainHead,iChainTail);
            }
            else
            {
                confs.setHead(iChainTail,timeRange[1] );
                confs.join(iChainHead,tailChain.next);
                
            }
            
        }
        else
        {
            // copy back old configurations
            confs.copyDataToBuffer(buffer,timeRanges[0],iChainHead);
            confs.copyDataToBuffer(buffer,timeRanges[1],timeRanges[0][1] - timeRanges[0][0]+1);

        }
    
    return accept;
    }
};


moveHead::moveHead(int maxAdvanceLength_) :
maxAdvanceLength(maxAdvanceLength_),buffer(maxAdvanceLength_,getDimensions()) , _levy(maxAdvanceLength)
{

}

bool moveHead::attemptRecedeMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG,int iChainHead,int l)
{

   auto headChain = confs.getChain(iChainHead);

  
   std::array<int,2> timeRange = {  headChain.head  - l, headChain.head-1};


   auto timeRanges=splitPeriodicTimeSlice(timeRange,confs.nBeads() );


   Real deltaS=0;

    auto & sPot= S.getPotentialAction();


    deltaS+=sPot.evaluate(confs,timeRanges[0],headChain.prev);
    deltaS+=sPot.evaluate(confs,timeRanges[1],iChainHead);

    bool accept = sampler.acceptLog(-deltaS,randG);

    if ( accept)
    {
        if ( timeRange[0] < 0 )
        {
            confs.setHead(headChain.prev,timeRanges[0][0]);
            confs.removeChain( iChainHead );
        }
        else
        {
            confs.setHead(headChain.prev,timeRange[0]);
        }
        
    }
    else
    {

    }


    return accept;
}


bool moveHead::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG)
{
    int iChainHead= std::floor(distr(randG) * confs.heads().size() );
    auto headChain = confs.getChain(iChainHead);
    int l = std::floor(distr(randG) * 2 *confs.nBeads() -confs.nBeads() );


    if (l>=0 )
    {
       return  attemptAdvanceMove(confs,S,randG,iChainHead,l);
    }
    else
    {
        throw attemptRecedeMove(confs,S,randG,iChainHead,std::abs(l) );
    }

    return false;
}



bool moveHead::attemptAdvanceMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG,int iChainHead,int l)
{

   auto headChain = confs.getChain(iChainHead);

  
   std::array<int,2> timeRange = { headChain.head-1, headChain.head + l};

   auto timeRanges = splitPeriodicTimeSlice( timeRange,confs.nBeads());

    Real deltaS=0;

    auto & sPot= S.getPotentialAction();

    // copy the head bead
    auto & data = confs.dataTensor();

    for(int d=0;d<getDimensions();d++)
    {
        tmpMean[d]=data(iChainHead,d,headChain.head -1);
    }

    // sample a new head 
    confSampler.sampleFreeParticlePosition(tmpPosition,tmpMean, S.getTimeStep() * l, randG  );

    for(int d=0;d<getDimensions();d++)
    {
        data(iChainHead,d,headChain.head + l)=tmpPosition[d];
    }

    // performs reconstruction

    _levy.apply(confs,timeRange,iChainHead,sPot.getTimeStep(),randG); 

    int iNewChain=-1;
    if ( timeRange[1] > confs.nBeads( ) )
    {
        std::array<int,2> newTimeRange={0,timeRange[1] - confs.nBeads() };

        int iNewChain=confs.pushChain(confs.getGroupByChain(iChainHead));
        confs.setHead(iChainHead,newTimeRange[1]);
        confs.setTail(iNewChain,-1);
        confs.copyData({confs.nBeads(),timeRange[1]},iChainHead,0,iNewChain);
        confs.setHead(iChainHead,confs.nBeads()  );
        deltaS+=sPot.evaluate(confs,newTimeRange,iNewChain);
    }
    else
    {
        confs.setHead(iChainHead,timeRange[1]);
        deltaS+=sPot.evaluate(confs,timeRange,iChainHead);
    }

    bool accept = sampler.acceptLog(-deltaS,randG);

    if (not accept)
    {
        if ( timeRange[1] > confs.nBeads() )
            {
                confs.removeChain(iNewChain);
            }
        confs.setHead(iChainHead,headChain.head);
    }
    else
    {
        if ( timeRange[1] > confs.nBeads() )
        {
            confs.join(iChainHead,iNewChain);
        }
    }

    confs.fillHead( iChainHead);
    confs.fillHead( iNewChain);

    return accept;
}


}