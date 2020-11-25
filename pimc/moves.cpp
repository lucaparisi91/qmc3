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

        auto & data = configurations.dataTensor() ;
        
        auto timeRanges = splitPeriodicTimeSlice(timeRange,configurations.nBeads() );

        const auto & currentChain = configurations.getChain(iChain);
        int iChainNext=currentChain.next;

        //int lastBeadChain= timeRange[1] > configurations.nBeads() ? iChainNext : iChain ;
        //int lastBeadTime= timeRange[1] > configurations.nBeads() ? timeRanges[1][1] : timeRange[1] ;


 /*         if (lastBeadChain == -1 )
        {
            throw invalidState("Crossing a head in reconstruction.");
        }

        for (int d=0;d<getDimensions();d++)
            {
                data(iChain,d,timeRange[1])=data(lastBeadChain,d,lastBeadTime);
            }
  */

        
        // performs the actual copy
        for (int t=0;t<l-1;t++)
        {
           
            

            for (int d=0;d<getDimensions();d++)
                {
                Real eta = gauss(randG);
                mean[d] = 
                (
                data(iChain,d,t+timeRange[0])*(l-t-1) 
                + data(iChain,d , l + timeRange[0]  ) )
                /( l - t) 
                ;

                Real variance = (l-t-1) * 1. /(l-t) *timeStep;

                data(iChain,d, timeRange[0] + t + 1 ) = mean[d] + eta *sqrt(variance) ;
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



      if (timeRange[1] > confs.nBeads() )//copy the end bead in second chain to the first chain
    {
        confs.copyData( { timeRanges[1][1] , timeRanges[1][1] } , iChainNext, timeRange[1],iChain );
    }
    _levy.apply(confs,timeRange,iChain,S.getTimeStep(),randG);
    confs.copyData( { timeRanges[0][1]+1 , timeRange[1] } , iChain, 0,iChainNext ); // time periodic boundary conditions



    auto  sNew= S.evaluate(confs,timeRanges[0], iChain) ;
    sNew+=S.evaluate(confs,timeRanges[1], iChainNext) ;

    const auto actionDifference = sNew - sOld;

    bool accepted = sampler.acceptLog(-actionDifference,randG);

    if (! accepted)
    {
        // copy back old beads
        copyFromBuffer(confs,timeRanges[0],iChain);
        copyFromBuffer(confs,timeRanges[1],iChainNext, timeRanges[0][1]  - timeRanges[0][0] + 1 );
         confs.fillHead(iChain);
   
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

std::ostream & tableMoves::operator>> (std::ostream & os)
{
    for (int i=0;i<_moves.size();i++)
    {
        if (_nTrials[i] > 0 )
        {
            os << _names[i] << ":\t" << acceptanceRatio(i) << std::endl;
        }
        else
        {
            os << _names[i] << ":\t" << "-" << std::endl;
        }        
    }

    return os;

}

bool tableMoves::attemptMove(configurations_t & confs, firstOrderAction & S,randomGenerator_t & randG)
{
    int iMove=sample(randG);

    auto  move = (_moves[ iMove ]);
    bool success=move->attemptMove(confs,S,randG);

    _nTrials[iMove]++;

    if (success)
    {
    _nSuccess[iMove]++;
    }


    return success;
};


void tableMoves::push_back(move * move_,Real weight,const std::string & name)
{
    totalWeight+=weight;
    accumulatedWeights.push_back(totalWeight);
    _moves.push_back(move_);
    _names.push_back(name);
    _nTrials.push_back(0);
    _nSuccess.push_back(0);

};


int tableMoves::sample(randomGenerator_t & randG)
{
    int iMove = sampler.sample(accumulatedWeights,totalWeight,randG);
    return iMove;
};

openMove::openMove(Real C_ , int maxReconstructedLength_) : C(C_), _levy(maxReconstructedLength_) ,  _maxReconstructedLength(maxReconstructedLength_) ,buffer(maxReconstructedLength_,getDimensions()){}


translateMove::translateMove(Real max_delta, int maxBeads) : _max_delta(max_delta),buffer(maxBeads,getDimensions())  , distr(-1.,1.)
{


}

bool translateMove::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{

    auto & Spot = S.getPotentialAction();
    std::array<int,2> timeRange={0,confs.nBeads()-1};

    // sample a chain
    int iChain = confSampler.sampleChain(confs,randG);
    auto & data = confs.dataTensor();

    auto chainsInThePolimer = confs.buildPolimerList(iChain); 

    Real deltaS = 0;

    // evaluate the action in the old configuration

     for (auto iCurrentChain : chainsInThePolimer)
     {
        deltaS-= Spot.evaluate(confs,timeRange,iCurrentChain);   
     }




    int iSeq=0; // ith chain in the list
    for (auto iCurrentChain : chainsInThePolimer)
    {

        // save old data to buffer
        confs.copyDataToBuffer(buffer,{0,confs.nBeads()},iCurrentChain,iSeq*confs.nBeads());


         // translate the whole 
        for (int d=0;d<getDimensions();d++)
        {
            delta[d]=distr(randG)*_max_delta;
        }

        for(int t=0;t<=confs.nBeads();t++)
            {
                for (int d=0;d<getDimensions();d++)
                {
                    data(iCurrentChain,d,t)+=delta[d];
                }    

            }

        iSeq++; 
    }


    //evaluate the action in the new configuration

     for (auto iCurrentChain : chainsInThePolimer)
     {
        deltaS+= Spot.evaluate(confs,timeRange,iCurrentChain);   
     }



    bool accept = sampler.acceptLog(-deltaS,randG);

    if ( not accept)
    {
        // copy back configurations in al chains of the permutation cycle
        int iSeq=0;
        for (auto iCurrentChain : chainsInThePolimer)
        {
          confs.copyDataFromBuffer(buffer,{0,confs.nBeads()},iCurrentChain,iSeq*confs.nBeads());
        }
        iSeq++;
    }

    return accept;

}

bool openMove::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{

    
    int iChain = confsSampler.sampleChain(confs,randG);
    int iChainNext=confs.getChain(iChain).next;

    Real var=2*D*timeStep;

    int iHead= std::floor( uniformRealNumber(randG) * ( confs.nBeads() ) );
    
    int l= std::floor( uniformRealNumber(randG) * (_maxReconstructedLength -2) ) + 1 ; // distance from itime where the head is formed

    
    // if time range is negative shift to previous chain and translate the time slice by the number of beads
    std::array<int,2> timeRange{iHead - l - 1,iHead};
    if ( timeRange[0] <  0 )
    {
        timeRange[0]+=confs.nBeads();
        timeRange[1]+=confs.nBeads();
        iChain=confs.getChain(iChain).prev;
    }
    int iChainHead=  timeRange[1] > confs.nBeads() ? iChainNext : iChain;
    auto deltaS=0;

    auto timeRanges = splitPeriodicTimeSlice(timeRange,confs.nBeads());
  
    auto & data = confs.dataTensor();
    

    std::array<Real,3> difference;
    std::array<Real,3> oldBead;

    deltaS-=S.evaluate(confs,timeRanges[0],iChain);
    deltaS-=S.evaluate(confs,timeRanges[1],iChainNext);


    confs.copyDataToBuffer(buffer,timeRanges[0],iChain,0);
    confs.copyDataToBuffer(buffer,timeRanges[1],iChainNext,timeRanges[0][1]- timeRanges[0][0] + 1);

    for(int d=0;d<getDimensions();d++)
    {
        oldBead[d]=data(iChainHead,d,iHead);
    }


    // generates the head
    Real distanceSquared=0;
    const auto & geo = S.getGeometry();
    std::array<Real,3> headPosition;
    std::array<Real,3> startPosition;

    for (int d=0;d<getDimensions();d++)
    {
        startPosition[d]=data(iChain,d,timeRange[0]);
    }

    Real mass = confs.getGroupByChain(iChain).mass;

    confsSampler.sampleFreeParticlePosition(headPosition,startPosition,timeStep*(l+1),randG,mass);

    for (int d=0;d<getDimensions();d++)
    {
        data(iChain,d,timeRange[1])=headPosition[d];
    }

    
    // perform levy reconstruction on l beads
    _levy.apply(confs,timeRange,iChain,timeStep,randG);
    confs.copyData({confs.nBeads(),timeRange[1]},iChain,0,iChainNext);

    // evaluates the action


    deltaS+=S.evaluate(confs,timeRanges[0],iChain);
    deltaS+=S.evaluate(confs,timeRanges[1],iChainNext);


    // compute the acceptance ratio
    for (int d=0;d<getDimensions();d++)
    {
    
        difference[d]=
            geo.difference( 
                data(iChain,d,timeRange[1])-data(iChain,d,timeRange[0]),d
            );
    }
    auto propRatio = -deltaS - freeParticleLogProbability(difference,S.getTimeStep(),mass);

    bool accept = sampler.acceptLog(propRatio,randG);

    if ( accept)
    {

        if (iHead==confs.nBeads() ) // next chain will be the new tail
        {
            confs.setHead(iChain,iHead);
        }
        else // creates a new chain to contain the new tail
        {
            
            int tailChain=confs.pushChain(confs.getGroupByChain(iChain));
            confs.setTail(tailChain,iHead - 1);
            confs.setHead(tailChain,confs.nBeads() );
            confs.setHead(iChainHead,iHead);
            confs.join(tailChain,confs.getChain(iChainHead).next);
            // copy the upper chain in the new tail, including the old value of the head
            confs.copyData({iHead+1,confs.nBeads()-1}, iChainHead, tailChain );
            for (int d=0;d<getDimensions();d++)
            {
                data(tailChain,d,iHead)=oldBead[d];
            }
            confs.fillHead(tailChain);
        }
    }
    else
    {

    confs.copyDataFromBuffer(buffer,timeRanges[0],iChain,0);
    confs.copyDataFromBuffer(buffer,timeRanges[1],iChainNext,timeRanges[0][1]- timeRanges[0][0] + 1);
    confs.fillHead(iChain);
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