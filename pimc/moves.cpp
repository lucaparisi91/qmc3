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

    void levyReconstructor::apply (configurations_t & configurations, std::array<int,2> timeRange,int iChain ,const action & S,randomGenerator_t & randG)
    {
        /* Does not take into account time periodic boundary conditions
        Reconstruct between timeRange[1] -1 and timeRange[0] + 1
        timeRange[0] and timeRange[1] are used as boundary conditions
        */

        Real timeStep=S.getTimeStep();
        const auto & geo = S.getGeometry();

        int l = timeRange[1] - timeRange[0];

        auto & data = configurations.dataTensor() ;


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

        std::array<Real,getDimensions()> xBottom;
        std::array<Real,getDimensions()> xTop;
        
        // performs the actual copy
        for (int t=0;t<l-1;t++)
        {

            for (int d=0;d<getDimensions();d++)
                {

                Real eta = gauss(randG);

                //xBottom[d]=geo.difference( data(iChain,d,t+timeRange[0]) - data(iChain,d,t+1+timeRange[0]),d) + data(iChain,d,t+1+timeRange[0]) ; 
                
                //xTop[d]=geo.difference( data(iChain,d,l+timeRange[0]) - data(iChain,d,t+1+timeRange[0]),d) + data(iChain,d,t+1+timeRange[0]) ; 

                xBottom[d]=data(iChain,d,t+timeRange[0]);
                xTop[d]=data(iChain,d,l+timeRange[0]);






                mean[d] = 
                (
                xBottom[d]*(l-t-1) 
                + xTop[d] )
                /( l - t) 
                ;

                Real variance = (l-t-1) * 1. /(l-t) *timeStep;

                data(iChain,d, timeRange[0] + t + 1 ) = mean[d] + eta *sqrt(variance) ;
                }
        }
        

        //configurationsNew.copyData(  {configurationsNew.nBeads() , timeRange[1]-1} , iChain , 0  , iChainNext);
        //configurationsNew.fillHead(iChain);

    }

levyMove::levyMove(int maxBeadLength_) : _levy(maxBeadLength_) , uniformRealNumber(0,1),maxBeadLength(maxBeadLength_) , buffer(( maxBeadLength_+1)*2,getDimensions() ) {}

bool levyMove::attemptMove( configurations_t & confs, firstOrderAction & ST,randomGenerator_t & randG)
{

    int nChains = confs.nChains();
    int nBeads = confs.nBeads();
    auto & S = ST.getPotentialAction();

    auto & geo = S.getGeometry();


    int iChain=confsSampler.sampleChain(confs,randG);

    auto timeRange = tGen(randG,nBeads,maxBeadLength);

    auto timeRanges = splitPeriodicTimeSlice(timeRange,confs.nBeads());
    const auto & currentChain = confs.getChain(iChain);

    int iChainNext = currentChain.next;

    if (
        (timeRange[1] > confs.nBeads() and  (currentChain.next == -1 ) )
        or (   ( currentChain.next != -1 ) and
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
     confs.copyDataToBuffer(buffer,timeRanges[0],iChain);
     if (iChainNext!=-1)
     {
        confs.copyDataToBuffer(buffer,timeRanges[1],iChainNext, timeRanges[0][1]  - timeRanges[0][0] + 1);
     }

      if (timeRange[1] > confs.nBeads() )//copy the end bead in second chain to the first chain
    {
        confs.copyData( { timeRanges[1][1] , timeRanges[1][1] } , iChainNext, timeRange[1],iChain );

        for(int d=0;d<getDimensions();d++)
        {
            data(iChain,d,timeRange[1])-= data(iChainNext,d,0) - data(iChain,d,confs.nBeads() );
        }// ensures to reconstruct along a continuos path

        
    }
    _levy.apply(confs,timeRange,iChain,S,randG);
    if (iChainNext!=-1)
     {
         confs.copyData( { timeRanges[0][1]+1 , timeRange[1]-1 } , iChain, 0,iChainNext ); // time periodic boundary conditions

        for(int t=0;t<timeRanges[1][1];t++)
            for(int d=0;d<getDimensions();d++)
            {
                data(iChainNext,d,t)+=
                data(iChainNext,d,timeRanges[1][1]) - data(iChain,d,timeRange[1] );


                //std::cout << data(iChainNext,d,timeRanges[1][1]) - data(iChain,d,timeRange[1] )<<std::endl;
            }// ensures to reconstruct along a continuos path
     }




    auto  sNew= S.evaluate(confs,timeRanges[0], iChain) ;
    sNew+=S.evaluate(confs,timeRanges[1], iChainNext) ;

    const auto actionDifference = sNew - sOld;

    bool accepted = sampler.acceptLog(-actionDifference,randG);

    if (! accepted)
    {
        // copy back old beads
        confs.copyDataFromBuffer(buffer,timeRanges[0],iChain);
        if (iChainNext!=-1)
         {
        confs.copyDataFromBuffer(buffer,timeRanges[1],iChainNext, timeRanges[0][1]  - timeRanges[0][0] + 1);
         }
         
    }

    
    return accepted;
}

swapMove::swapMove(int maxStepLength_,int maxN) :
 maxStepLength(maxStepLength_),buffer(maxStepLength_*2,getDimensions() )
 , uniformRealNumber(0,1) ,_levy(maxStepLength_+2),particleSampler(maxN)
{}

bool swapMove::attemptMove(configurations_t & confs, firstOrderAction & S,randomGenerator_t & randG)
{
    if ( ! confs.isOpen() )
    {
        throw invalidState("Swap move can only be done in the G sector.");
    }
     auto & data = confs.dataTensor();
    // selects an head at random

    int _iChainHead =  std::floor(  uniformRealNumber(randG) * confs.heads().size() ) ;
    int iChainHead =  confs.heads()[_iChainHead] ;

    auto & Spot = S.getPotentialAction();

    // selects a time slice length at random. Refuse if reconstructed time slice crosses over the end bead

    int l = std::floor(uniformRealNumber(randG) * maxStepLength) + 1;
    //int l = maxStepLength;

    
    int iHead=confs.nBeads();

    
    // tower sampling a particle i with gaussian weights on the relative distances
   
    const auto & group = confs.getGroup( iChainHead );
    auto & geo = S.getGeometry();

    std::array<Real, 3> distance;
    particleSampler.reset();

    Real weightForwardMove = 0;
    
    for(int i=group.iStart;i<=group.iEnd;i++)
    {
        
        for(int d=0;d<getDimensions();d++)
        {
            distance[d]=geo.difference(  data(iChainHead,d,iHead) - data(i,d,l),d);
        }

        particleSelectionWeight=exp(freeParticleLogProbability(distance,S.getTimeStep()*l,group.mass));
        weightForwardMove+=particleSelectionWeight;
        particleSampler.accumulateWeight(particleSelectionWeight);
    }

    int iPartner=particleSampler.sample(randG) + group.iStart;

    if ( confs.getChain(iPartner).hasTail() )
    {
        return false;
    }

    // metropolic test based on ratio of forward and backward move
    Real weightBackwardMove= 0;

    for(int i=group.iStart;i<=group.iEnd;i++)
    {
        Real norm=0;
        for(int d=0;d<getDimensions();d++)
        {
            distance[d]=geo.difference(  data(iPartner,d,0) - data(i,d,l),d);
        }


        weightBackwardMove+=exp(freeParticleLogProbability(distance,S.getTimeStep()*l,group.mass));
    }


    

    Real deltaS=-(log(weightForwardMove) - log(weightBackwardMove)) ;

    
    const auto  partnerChain = confs.getChain(iPartner);

    
    deltaS-=Spot.evaluate(confs,{0,l-1},iPartner);

    confs.copyDataToBuffer(buffer,{0,l-1},iPartner);


    // performs levy reconstruction between the head and the bead
    for(int d=0;d<getDimensions();d++)
        {
            data(iPartner,d,0)=//data(iChainHead,d,iHead);
            geo.difference (- data (iPartner,d,l) + data(iChainHead,d,iHead) , d);
            data(iPartner,d,0)+=data (iPartner,d,l);
        }

    _levy.apply(confs,  {0,l },  iPartner, S ,randG );
    deltaS+=Spot.evaluate(confs,{0,l-1},iPartner);
    
    bool accept = metropolisSampler.acceptLog(-deltaS,randG);
    if (accept) 
    {
        confs.join(iChainHead,iPartner );
    }
    else
    {
        confs.copyDataFromBuffer(buffer,{0,l-1},iPartner);
    }    

    return accept;
}

std::ostream & sectorTableMoves::operator>> (std::ostream & os)
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

bool sectorTableMoves::attemptMove(configurations_t & confs, firstOrderAction & S,randomGenerator_t & randG)
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


void sectorTableMoves::push_back(move * move_,Real weight,const std::string & name)
{
    sampler.accumulateWeight(weight);
    _moves.push_back(move_);
    _names.push_back(name);
    _nTrials.push_back(0);
    _nSuccess.push_back(0);

};

int sectorTableMoves::sample(randomGenerator_t & randG)
{
    int iMove = sampler.sample(randG);
    return iMove;
};

openMove::openMove(Real C_ , int maxReconstructedLength_) : C(C_), _levy(maxReconstructedLength_+2) ,  _maxReconstructedLength(maxReconstructedLength_+2) ,buffer(2*(maxReconstructedLength_+2),getDimensions()),
gauss(0,1),uniformRealNumber(0,1){}

translateMove::translateMove(Real max_delta, int maxBeads) : _max_delta(max_delta),buffer(maxBeads+1,getDimensions())  , distr(-1.,1.)
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



    // translate the whole 
    for (int d=0;d<getDimensions();d++)
    {
        delta[d]=(distr(randG))*_max_delta;
    }

    int iSeq=0; // ith chain in the list
    for (auto iCurrentChain : chainsInThePolimer)
    {

        // save old data to buffer
        confs.copyDataToBuffer(buffer,{0,confs.nBeads()},iCurrentChain,iSeq*(confs.nBeads()+1));



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
          confs.copyDataFromBuffer(buffer,{0,confs.nBeads()},iCurrentChain,iSeq*(confs.nBeads()+1));
          iSeq++;
        }
    }

    return accept;

}

bool openMove::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{
    Real timeStep = S.getTimeStep();

    if ( confs.isOpen() )
    {
        throw invalidState("The configuration is already open.");
    }

    int iChain = confsSampler.sampleChain(confs,randG);
    int iChainTail=confs.getChain(iChain).next;

    int iHead=confs.nBeads();
    
    //int l= std::floor( uniformRealNumber(randG) * (_maxReconstructedLength -2) ) + 1 ; // distance from itime where the head is formed

    int l = _maxReconstructedLength - 2;

    int t0= iHead - l;
    Real deltaS=0;
    std::array timeRange={t0  , iHead -1 };

    auto & data = confs.dataTensor();

    std::array<Real,3> difference;
    auto & sPot = S.getPotentialAction();

    deltaS-=sPot.evaluate(confs,timeRange,iChain);

    confs.copyDataToBuffer(buffer,{t0,iHead},iChain,0);

    // generates the head
    
    const auto & geo = S.getGeometry();
    std::array<Real,getDimensions()> headPosition;
    std::array<Real,getDimensions()> startPosition;

    for (int d=0;d<getDimensions();d++)
    {
        startPosition[d]=data(iChain,d,t0);
    }


    for (int d=0;d<getDimensions();d++)
    {
        difference[d]=
               geo.difference( 
                data(iChainTail,d,0)-data(iChain,d,t0),d
            );

              //data(iChain,d,iHead)-data(iChain,d,t0);

        if (
            std::abs(data(iChain,d,iHead)-data(iChain,d,t0) ) > geo.getLBox(d)*0.5 )
            {
                std::cout <<    data(iChain,d,iHead)-data(iChain,d,t0) << std::endl;
                return false;
            }
        
        
   
    }


    
    Real mass = confs.getGroupByChain(iChain).mass;

    confsSampler.sampleFreeParticlePosition(headPosition,startPosition,timeStep*l,randG,mass);

    for (int d=0;d<getDimensions();d++)
    {
         if (
            std::abs(headPosition[d]-startPosition[d] ) > geo.getLBox(d)*0.5 )
            {
                std::cout <<    data(iChain,d,iHead)-data(iChain,d,t0) << std::endl;
                return false;
            }

    }

    for (int d=0;d<getDimensions();d++)
    {
        data(iChain,d,iHead)=headPosition[d];
    
    }


    


    // perform levy reconstruction on l beads
    _levy.apply(confs,{t0,iHead},iChain,S,randG);

    // evaluates the action
    deltaS+=sPot.evaluate(confs,timeRange,iChain);



    // compute the acceptance ratio
  
    auto propRatio = -deltaS - freeParticleLogProbability(difference,S.getTimeStep()*l,mass) + log(C);


    bool accept = sampler.acceptLog(propRatio,randG);

    if ( accept)
    {

        confs.setHead(iChain,iHead);
    }
    else
    {
        confs.copyDataFromBuffer(buffer,{t0,iHead},iChain,0);
    }

    return accept;
};


closeMove::closeMove(Real C_ , int maxReconstructionLength) : C(C_),_levy(2*(maxReconstructionLength+2)),_maxLength(maxReconstructionLength+2),buffer((maxReconstructionLength+2)*2,getDimensions()),gauss(0,1),uniformRealNumber(0,1)
{}

bool closeMove::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{

    int iChainHead=confs.heads()[std::floor(uniformRealNumber(randG) * confs.heads().size() )];
    int iChainTail=confs.tails()[std::floor(uniformRealNumber(randG) * confs.tails().size() )];
   
    auto timeStep = S.getTimeStep();

    auto & geo = S.getGeometry();

    int iHead=confs.nBeads();
    int iTail = 0;

    
    
    //int l= std::floor( uniformRealNumber(randG) * (_maxLength -2) ) + 1 ; // distance from itime where the head is formed
    int l = _maxLength - 2;
    
    auto & sPot = S.getPotentialAction();


    Real deltaS=0;
    int t0= iHead - l;

    std::array timeRange={iHead - l , iHead -1 };

    auto & data = confs.dataTensor();

    //std::cout << data(iChainHead,0,iHead)-data(iChainTail,0,iTail) << std::endl;


    std::array<Real,3> difference;

    deltaS-=sPot.evaluate(confs,{t0,iHead-1},iChainHead);

    confs.copyDataToBuffer(buffer,{t0,iHead },iChainHead,0);

    // copy first valid bead of tail of the head
    Real distanceSquared=0;
    std::array<Real,3> headPosition;
    std::array<Real,3> startPosition;
    std::array<Real,3> distance;


    for (int d=0;d<getDimensions();d++)
    {
        //std::cout << "close: " <<geo.difference( - data(iChainHead,d,t0) + data(iChainTail,d,iTail),d) << std::endl;


        data(iChainHead,d, iHead )=geo.difference( - data(iChainHead,d,t0) + data(iChainTail,d,iTail),d);
        //data(iChainHead,d, iHead )= - data(iChainHead,d,t0) + data(iChainTail,d,iTail);

        data(iChainHead,d, iHead )+=data(iChainHead,d,t0);

        
        //std::cout<<  ( - data(iChainHead,d,t0) + data(iChainHead,d,iHead),d) << std::endl;

        //std::cout <<   - data(iChainHead,d,0) + data(iChainHead,d,iHead) << std::endl;






    }


    Real mass = confs.getGroupByChain(iChainHead).mass;
    

    // perform levy reconstruction on l beads
    _levy.apply(confs,{t0,iHead},iChainHead,S,randG);

    // evaluates the action
    deltaS+=sPot.evaluate(confs,{t0,iHead-1},iChainHead);

    // compute the acceptance ratio
    for (int d=0;d<getDimensions();d++)
    {
        difference[d]=
            //geo.difference( 
            //    data(iChainHead,d,t0)-data(iChainHead,d,iHead),d
            //);
         data(iChainHead,d,t0)-data(iChainHead,d,iHead);
    }

/* 
     int winding =  (data(iChainHead,0,iHead)-data(iChainTail,0,iTail))/2;
     
     {
        std::cout << winding << std::endl;
     } */



    auto propRatio = -deltaS + freeParticleLogProbability(difference,S.getTimeStep()*l,mass) -log(C);


    //std::cout << propRatio << std::endl;


    bool accept = sampler.acceptLog(propRatio,randG);

    if ( accept)
    {
        confs.join(iChainHead,iChainTail);
        //confs.fillHead(iChainHead);

        /* for(int d=0;d<getDimensions();d++)
        {
            std::cout<< - data(iChainHead,d,iHead) + data(iChainTail,d,iTail) << " ";
        }
        std::cout << std::endl; */
        
    }
    else
    {
        confs.copyDataFromBuffer(buffer,{t0,iHead },iChainHead,0);
    }
    

    

    return accept;

};


moveHead::moveHead(int maxAdvanceLength_) :
_maxReconstructedLength(maxAdvanceLength_+2),buffer((maxAdvanceLength_+2)*2,getDimensions()) , _levy((maxAdvanceLength_+2)*2),gauss(0,1),uniformRealNumber(0,1)
{

}


moveTail::moveTail(int maxAdvanceLength_) :
_maxReconstructedLength(maxAdvanceLength_+2),buffer((maxAdvanceLength_+2)*2,getDimensions()) , _levy((maxAdvanceLength_+2)*2),gauss(0,1),uniformRealNumber(0,1)
{
    
}


bool moveHead::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG)
{
    Real timeStep = S.getTimeStep();


    int iChain=confs.heads()[std::floor(uniformRealNumber(randG) * confs.heads().size() )];
    

    int iHead=confs.nBeads();


    
    int l= std::floor( uniformRealNumber(randG) * (_maxReconstructedLength -2) ) + 1 ; // distance from itime where the head is formed

    int t0= iHead - l;
    Real deltaS=0;
    std::array timeRange={t0  , iHead -1 };

    auto & data = confs.dataTensor();

    std::array<Real,3> difference;
    
    auto & sPot = S.getPotentialAction();

    deltaS-=sPot.evaluate(confs,timeRange,iChain);

    confs.copyDataToBuffer(buffer,{t0,iHead},iChain,0);

    // generates the head
    Real distanceSquared=0;
    const auto & geo = S.getGeometry();
    std::array<Real,getDimensions()> headPosition;
    std::array<Real,getDimensions()> startPosition;

    for (int d=0;d<getDimensions();d++)
    {
        startPosition[d]=data(iChain,d,t0);
    }

    Real mass = confs.getGroupByChain(iChain).mass;
    Real var=2*D*timeStep/ mass;

    confsSampler.sampleFreeParticlePosition(headPosition,startPosition,timeStep*l,randG,mass);

    for (int d=0;d<getDimensions();d++)
    {
        data(iChain,d,iHead)=headPosition[d];
    }
    
    // perform levy reconstruction on l beads
    _levy.apply(confs,{t0,iHead},iChain,S,randG);

    // evaluates the action
    deltaS+=sPot.evaluate(confs,timeRange,iChain);


    auto propRatio = -deltaS;

    bool accept = sampler.acceptLog(propRatio,randG);

    if ( accept)
    {

        
    }
    else
    {
        confs.copyDataFromBuffer(buffer,{t0,iHead},iChain,0);
    }

    return accept;

}

bool moveTail::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG)
{
    Real timeStep = S.getTimeStep();

    int iChain=confs.tails()[std::floor(uniformRealNumber(randG) * confs.tails().size() )];

    if (not confs.getChain(iChain).hasTail() )
    {
        throw invalidState("No tail to move.");
    }

    int iTail=0;

    
    int l= std::floor( uniformRealNumber(randG) * (_maxReconstructedLength -2) )  + 1; // distance from itime where the head is formed
    //int l=_maxReconstructedLength - 2;


    int t1= iTail + l;
    Real deltaS=0;
    std::array timeRange={iTail  , t1 - 1 };


    auto & data = confs.dataTensor();

    std::array<Real,3> difference;
    
    auto & sPot = S.getPotentialAction();

    deltaS-=sPot.evaluate(confs,timeRange,iChain);

    confs.copyDataToBuffer(buffer,{iTail,t1},iChain,0);

    // generates the tail
    Real distanceSquared=0;
    const auto & geo = S.getGeometry();
    std::array<Real,getDimensions()> tailPosition;
    std::array<Real,getDimensions()> startPosition;

    for (int d=0;d<getDimensions();d++)
    {
        startPosition[d]=data(iChain,d,t1);
    }

    Real mass = confs.getGroupByChain(iChain).mass;
    Real var=2*D*timeStep/ mass;

    confsSampler.sampleFreeParticlePosition(tailPosition,startPosition,timeStep*l,randG,mass);

    for (int d=0;d<getDimensions();d++)
    {
        data(iChain,d,iTail)=tailPosition[d];
    }

    
    // perform levy reconstruction on l beads
    _levy.apply(confs,{iTail,t1},iChain,S,randG);

    // evaluates the action
    deltaS+=sPot.evaluate(confs,timeRange,iChain);

    auto propRatio = -deltaS;

    bool accept = sampler.acceptLog(propRatio,randG);

    if ( accept)
    {

    }
    else
    {
        confs.copyDataFromBuffer(buffer,{iTail,t1},iChain,0);
    }

    return accept;

}

void tableMoves::push_back( move * move_,Real weight,sector_t sector,const std::string & name)
{
    if (sector == sector_t::diagonal)
    {
        closedTab.push_back(move_,weight,name);
    }
    else if (sector == sector_t::offDiagonal)
    {
        openTab.push_back(move_,weight,name);
    }

}


bool tableMoves::attemptMove(configurations_t & confs, firstOrderAction & S,randomGenerator_t & randG)
{
    if (  confs.isOpen() )
    {
        nOpenSectorMoves++;
        return openTab.attemptMove(confs,S,randG);
    }
    else
    {
        nClosedSectorMoves++;
        return closedTab.attemptMove(confs,S,randG);
    }
        
}

std::ostream & tableMoves::operator>> (std::ostream & os)
{
    os << "----------Open Sector" << "-----------" << std::endl;
    openTab >> os;
    os << "----------Closed Sector" << "-----------" << std::endl;
    closedTab >> os;
    os << "----------------" << std::endl;

    if ((nOpenSectorMoves + nClosedSectorMoves) > 0 )
    {
        os << "Open Sector fraction: " << nOpenSectorMoves/(nOpenSectorMoves + nClosedSectorMoves) << std::endl;
    }

    return os;
}


}


