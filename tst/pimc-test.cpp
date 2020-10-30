#include "gtest/gtest.h"
#include "geometryPMC.h"
#include "traits.h"
#include "tools.h"
#include "../pimc/action.h"
#include "../pimc/pimcConfigurations.h"
#include "../pimc/moves.h"
#include "../pimc/pimcObservables.h"


TEST(distances,updateSingleParticleDistances)
{
    pimc::geometryPBC_PIMC geo(10,10,10);
    int D = getDimensions();
    int N= 50;

    state_t positions(N,D);
    positions.setRandom();
    positions=positions - 0.5 ;

    auto difference = geo.differencesTwoBody(positions);
    positions(0,0)+=1;
    auto difference2 = geo.differencesTwoBody(positions);

    geo.differencesTwoBody(difference, positions , 0 );

    for (int k=0;k<difference2.rows();k++)
    {
        ASSERT_NEAR( difference2(k,0) , difference(k,0) , 1e-6  ) ;
    }

    int T = 50;
    

    Eigen::Tensor<Real, 3> data(N,getDimensions(),T);
    data.setRandom();

    Eigen::Tensor<Real,3> springDifferences(T,getDimensions() , N );


    geo.updateSpringDifferences(  springDifferences, data, {0,T-1} , {0,N-1} );

    for (int t=0;t<T;t++)
    {
        for (int i=0;i<N;i++)
        {
            int d=0;
            ASSERT_NEAR( springDifferences(t,d,i) , data(i,d,(t+1)%T) - data(i,d,t) , 1e-4 );
        }
    }

    Eigen::Tensor<Real,3> potentialDifferences(N ,getDimensions() , T );

    int jChain = N-1;
    geo.updateEqualTimeDifferences(potentialDifferences,data, {0,T-1}, jChain , {0,N-1});


    for (int t=0;t<T;t++)
    {
        for (int i=0;i<N;i++)
        {
            for(int d=0;d<getDimensions();d++)
            {
            ASSERT_NEAR( potentialDifferences(i,d,t) , data(i,d,t) - data(jChain,d,t) ,1e-4);
            }
        }
    }


    //auto springDistances = geo.springDistances( timeConfigurations );

    //geo.springDistances(springDistances,timeConfigurations,0,2 , 0, 3  );

}


TEST(configurations, init)
{
    const int N = 20;
    const int M = 30;

    pimc::particleGroup groupA{ 0 , N-1, N , 1.0};
    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});

    auto & data = configurations.dataTensor();

    data.setRandom();

    pimc::geometryPBC_PIMC geo(10,10,10);

    Real timeStep = 1e-2;



    pimc::kineticAction sT(timeStep, configurations.nChains(), configurations.nBeads()  , geo);



    

    Real currentKineticAction = sT.evaluate(configurations);

    Real kineticActionSimple=0;

    for(int n=0;n<N;n++)
    {
        for (int t=0;t<M;t++)
            {
                for(int d=0;d<getDimensions();d++)
                {
                    kineticActionSimple+=std::pow(data(n,d,(t+1)%M ) - data(n,d,t) ,2) /(4*0.5*timeStep) ;   
                }
            }
    }

    ASSERT_NEAR(currentKineticAction,kineticActionSimple,1e-5);

     auto harmonicPotential = [](Real x,Real y , Real z) {return x*x + y*y + z*z;};

     pimc::potentialActionOneBody<decltype(harmonicPotential)> pot1(timeStep,harmonicPotential,geo);

     auto v = pot1.evaluate(configurations);
     Real vCheck=0;

     for (int t=0;t<M;t++)
        for(int n=0;n<N;n++)
        {
            for(int d=0;d<getDimensions();d++)
            {
                vCheck+=data(n,d,t) * data(n,d,t) ;
            }
        }

        vCheck*=timeStep;
    ASSERT_NEAR(v,vCheck,1e-4);

    pimc::potentialActionTwoBody<decltype(harmonicPotential)> pot2(timeStep, N , M , harmonicPotential,geo) ;

    int jChain = 3;

    v=pot2.evaluate(configurations, {0,M-1} , jChain ,  {0,N-1} );
    vCheck=0;

    for ( int t=0;t<M;t++)
        for ( int i=0; i<N;i++ )
        {
            if (jChain != i)
            {
            for (int d=0;d<getDimensions();d++)
            {
                vCheck+= std::pow( data(i,d,t) - data(jChain,d,t) ,2);
            }
            }
        }

    ASSERT_NEAR(v,vCheck,1e-3);

    v=pot2.evaluate(configurations);

    vCheck=0;
    for ( int t=0;t<M;t++)
        for ( int i=0; i<N;i++ )
            for(int j=0;j<i;j++)
            {
                for(int d=0;d<getDimensions();d++)
                {
                    vCheck+= std::pow( data(i,d,t) - data(j,d,t) ,2);
                }

            }
    ASSERT_NEAR(v,vCheck,1e-3);



}


TEST(action , evaluation)
{
    const int N = 100;
    const int M = 30;

    pimc::particleGroup groupA{ 0 , N-1, N,  1.0};

    pimc::pimcConfigurations configurations(M, getDimensions() , {groupA});

    auto & data = configurations.dataTensor();

    data.setRandom();

    pimc::geometryPBC_PIMC geo(10,10,10);


    Real timeStep = 1e-2;

    pimc::kineticAction sT(timeStep, configurations.nChains(), configurations.nBeads(), geo);

    

    auto oldKineticAction = sT.evaluate(configurations);

    randomGenerator_t randG;

    std::array< std::array<int,2> ,2 > timeSlices ;
    timeSlices[0]={25,29};
    timeSlices[1]={0,9};

    int iChain = 55;




    auto sOld = sT.evaluate(configurations,timeSlices,iChain);

    auto deltaSOld= sT.evaluate(configurations, timeSlices, iChain);



    std::uniform_real_distribution<Real> unifDis(0,1);
    Real delta=0.5;

    for (int t=timeSlices[0][0];t<=timeSlices[0][1];t++)
    {
        for(int d=0;d<=getDimensions();d++)
        {
            data(iChain,d,t) += delta * unifDis(randG); 
        }
        
    }


    auto deltaSNew = sT.evaluate(configurations,timeSlices,iChain);

    auto sNew = sT.evaluate(configurations,timeSlices,iChain);


    auto deltaS = deltaSNew - deltaSOld;
    auto deltaSCheck = sNew - sOld;

    ASSERT_NEAR(deltaS,deltaSCheck,1e-4);

}


TEST(moves,levy_reconstructor)
{
    int seed = 30;
    const int N = 50;
    const int M = 30;

    pimc::particleGroup groupA{ 0 , N-1, N ,  1.0};
    randomGenerator_t randG(seed+5);

    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});

    pimc::pimcConfigurations configurations2(M , getDimensions() , {groupA});

    Real timeStep=1e-3;

    std::array<int,2> timeSlice= {10,26};


    pimc::levyReconstructor levy(timeStep);
    int iChain = 30;

    configurations.dataTensor().setRandom();
    configurations2=configurations;


    levy.apply(configurations,configurations,iChain,timeSlice,randG);

    auto &  data = configurations.dataTensor();
    auto &  data2 = configurations2.dataTensor();

     for (int i=0;i<N;i++)
    {
        for (int t=0; t<=timeSlice[0];t++)
         {
            for (int d=0;d<getDimensions();d++)
            ASSERT_NEAR( data2(i,d,t) , data(i,d,t) ,1e-4);
        }

        for (int t=timeSlice[1]; t<M;t++)
         {
            for (int d=0;d<getDimensions();d++)
            ASSERT_NEAR( data2(i,d,t) , data(i,d,t) ,1e-4);
        }

        if (i == iChain)
        {
            for (int t=timeSlice[0]+1; t<timeSlice[1];t++)
         {
            for (int d=0;d<getDimensions();d++)
            ASSERT_FALSE( abs( data2(i,d,t) - data(i,d,t) ) < 1e-5);
        }

        };



    }
 

    // test if gaussian distribution is sampled for a 3 point time slice

    int iSampled=29;
    iChain=3;
    timeSlice={iSampled-1,iSampled+1};
    Real alpha=1./(4*0.5*M_PI*timeStep);

    std::array<Real,getDimensions()> averagePosition {0,0,0};
    averagePosition={0,0,0};

    std::array<Real,getDimensions()> averagePositionSquared {0,0,0};
    std::array<Real,getDimensions()> averagePositionExtremes {0,0,0};
    int nSteps=10000;

    for( int d=0;d<getDimensions();d++)
    {
        averagePositionExtremes[d]=0.5*( data(iChain,d,iSampled-1) + data(iChain,d,(iSampled+1)%M) );
    }

    for (int i=0;i< nSteps ; i++)
    {
        levy.apply(configurations,configurations,iChain,timeSlice,randG);

        for(int d=0;d<getDimensions();d++)
        {
            averagePosition[d]+=data(iChain,d,iSampled);
            averagePositionSquared[d]+=data(iChain,d,iSampled) * data(iChain,d,iSampled) ;
        }
    }


    for (int d=0;d<getDimensions();d++)
    {
        averagePosition[d]/=nSteps;
        averagePositionSquared[d]/=nSteps;

        Real sigma= averagePositionSquared[d] - averagePosition[d]*averagePosition[d];
        Real sigmaExpected = timeStep/2.;

        Real meanDeviation=std::abs((averagePosition[d]  -  averagePositionExtremes[d])/averagePositionExtremes[d]);
        Real sigmaDeviation =  std::abs(sigmaExpected - sigma)/sigmaExpected;

        ASSERT_LT(sigmaDeviation,1e-2);
        ASSERT_LT(meanDeviation,1e-3);

    }



}

TEST(moves,levy)
{
    int seed = 30;
    const int N = 10;
    const int M = 50;

    pimc::particleGroup groupA{ 0 , N-1, N, 1.0};
    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});
     auto & data = configurations.dataTensor();
    data.setRandom();

    pimc::geometryPBC_PIMC geo(10,10,10);

    randomGenerator_t randG(seed);

    pimc::pimcConfigurations configurations2(M , getDimensions() , {groupA});

    
    Real timeStep=1e-3;
    std::array<int,2> timeSlice= {10,26};

    pimc::kineticAction sT(timeStep, N + 1 , M  , geo);

    auto V = [](Real x, Real y , Real z) {return 0.5*(x*x + y*y + z*z );};
    

    pimc::potentialActionOneBody<decltype(V)> sV(timeStep,V ,geo);

    pimc::firstOrderAction S(&sT, & sV);

    pimc::levyReconstructor levy(timeStep);
    pimc::levyMove mover(levy,20);
    int success = 0;

    auto & data2 = configurations2.dataTensor();

    int nSteps=10000;
    for (int i=0;i<nSteps;i++)
    {
        configurations2=configurations;
        bool accept = mover.attemptMove(configurations,S,randG);

        if (accept)
        {
            success+=1;
        }
        else
        {
            for(int t=0;t<M;t++)
                for(int i=0;i<N;i++)
                {
                    for(int d=0;d<getDimensions();d++)
                    {
                        ASSERT_NEAR( data(i,d,t) , data2(i,d,t) ,1e-3);
                        ASSERT_GT( std::abs(data(i,d,t)) , 0 );
                    }
                }
        }

    }

    ASSERT_GT(success,0);
    EXPECT_LT(success,nSteps);

 }



TEST(configurations, io)
{

    const int N = 10;
    const int M = 50;

    pimc::particleGroup groupA{ 0 , N-1, N , 1.0};
    
    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});
    
    auto & data = configurations.dataTensor();
    data.setRandom();
    int time=10;

    int iWorm=configurations.open(time, 5 );

    ASSERT_EQ( configurations.worms().size() , 1 );
    ASSERT_EQ(iWorm,0);

    auto worm = configurations.worms()[iWorm];

    
    configurations.save("testConfig");

    pimc::pimcConfigurations configurations2;

    configurations2.load("testConfig");

    ASSERT_EQ(configurations.nChains() ,configurations2.nChains() );
    ASSERT_EQ(configurations.nBeads() ,configurations2.nBeads() );

    const auto & data2 = configurations2.dataTensor();
    const auto & mask2 = configurations2.getMask();
    const auto & mask = configurations.getMask();

    for (int t=0;t<time;t++)
    {
            for(int d=0;d<getDimensions();d++)
            {
                ASSERT_NEAR(  data(worm.iChainHead,d,t) , data2(worm.iChainTail,d,t) , 1e-3 );

                ASSERT_EQ(  mask(t,worm.iChainTail) , 0);
                ASSERT_EQ(  mask(t,worm.iChainHead) , 1);

            }
    }

    for (int t=time+1;t<configurations.nBeads();t++)
    {
            for(int d=0;d<getDimensions();d++)
            {
                ASSERT_NEAR(  data(worm.iChainTail,d,t) , data2(worm.iChainTail,d,t) , 1e-3 );

                
                ASSERT_EQ(  mask(t,worm.iChainTail) , 1);
                ASSERT_EQ(  mask(t,worm.iChainHead) , 0);
              
            }
        
    }

    configurations.close(iWorm);

    for (int t=0;t<configurations.nBeads();t++)
    {
            for(int d=0;d<getDimensions();d++)
            {
                ASSERT_NEAR(  data(worm.iChainTail,d,t) , data2(worm.iChainTail,d,t) , 1e-3 );

                ASSERT_EQ(  mask(t,worm.iChainTail) , 1);

            }
    }

    ASSERT_EQ(configurations.worms().size(),0);


}

TEST(run,free_harmonic_oscillator)
{   
    int N=10;
    int M=100;

    pimc::geometryPBC_PIMC geo(30,30,30);

    Real timeStep = 1e-1;
    pimc::particleGroup groupA{ 0 , N-1, N , 1.0};
    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});

    configurations.dataTensor().setRandom();

    pimc::levyReconstructor reconstructor(timeStep);

    pimc::levyMove freeMoves(reconstructor, 10);

    pimc::tableMoves table;

    table.push_back(& freeMoves,timeStep);

    randomGenerator_t randG(100);

    pimc::kineticAction sT(timeStep, configurations.nChains() , M  , geo);

    auto V = [](Real x, Real y , Real z) {return 0.5*(x*x + y*y + z*z) ;};

    pimc::potentialActionOneBody<decltype(V)> sV(timeStep,V ,geo);

    pimc::firstOrderAction S(&sT, & sV);
    int nTimes = 10000;
    int success = 0;
    int subSteps=100;

    pimc::thermodynamicEnergyEstimator energyEstimator;
    Real e=0;
    Real e2=0;

    for (int i=0;i< nTimes ; i++)
    {
        for (int j=0;j<subSteps;j++)
        {
            auto & move = table.sample(randG);
            bool accepted=move.attemptMove(configurations, S, randG);

            if (accepted)
            {success+=1;}

        }
        
        Real tmp=energyEstimator(configurations,S);
        e+=tmp;
        e2+=tmp*tmp;

        //std::cout << e << std::endl;

        std::cout << "Energy: " << e/( (i+1)) << std::endl;
        std::cout << "Acceptance ratio: " << success*1./((i+1)*subSteps) << std::endl;
    }

    ASSERT_TRUE( (success*1./nTimes )> 0);
    e/=nTimes;
    e2/=nTimes;


    //std::cout << e << " " << std::sqrt(e2 - e*e) << std::endl;

}