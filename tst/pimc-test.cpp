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
#if DIMENSIONS == 3
TEST(configurations, init)
{
    const int N = 20;
    const int M = 30;

    pimc::particleGroup groupA{ 0 , N-1, N , 1.0};
    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});

    auto & data = configurations.dataTensor();

    data.setRandom();

    configurations.fillHeads();



    pimc::geometryPBC_PIMC geo(10,10,10);

    Real timeStep = 1e-2;



    pimc::kineticAction sT(timeStep, configurations.nChains(), configurations.nBeads()  , geo);



    

    Real currentKineticAction = sT.evaluate(configurations);

    Real kineticActionSimple=0;

    for(int n=0;n<N;n++)
    {
        for (int d=0;d<getDimensions();d++)
        {
             ASSERT_NEAR(data(n,d,M), data(n,d,0)  , 1e-5);
        }
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

#endif

#if DIMENSIONS == 3


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


    pimc::levyReconstructor levy(M);
    int iChain = 30;

    configurations.dataTensor().setRandom();
    configurations.fillHeads();

    configurations2=configurations;


    levy.apply(configurations,timeSlice,iChain,timeStep,randG);

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
    int nSteps=100000;

    for( int d=0;d<getDimensions();d++)
    {
        averagePositionExtremes[d]=0.5*( data(iChain,d,iSampled-1) + data(iChain,d,(iSampled+1)%M) );
    }

    for (int i=0;i< nSteps ; i++)
    {
        levy.apply(configurations,timeSlice,iChain,timeStep,randG);


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

#endif

TEST(moves,levy)
{
    int seed = 30;
    const int N = 10;
    const int M = 50;

    pimc::particleGroup groupA{ 0 , N-1, N, 1.0};
    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});
     auto & data = configurations.dataTensor();
    data.setRandom();
    configurations.fillHeads();

    pimc::geometryPBC_PIMC geo(10,10,10);

    randomGenerator_t randG(seed);

    pimc::pimcConfigurations configurations2(M , getDimensions() , {groupA});

    
    Real timeStep=1e-2;
    std::array<int,2> timeSlice= {10,26};

    pimc::kineticAction sT(timeStep, N + 1 , M  , geo);
  
    #if DIMENSIONS == 1
    auto V = pimc::makePotentialFunctor(
         [](Real x) {return 0.5*(x*x ) ;} ,
         [](Real x) {return 0.5*x ;} 
         );    
    #endif


    pimc::potentialActionOneBody<decltype(V)> sV(timeStep,V ,geo);

    pimc::firstOrderAction S(&sT, & sV);

    pimc::levyReconstructor levy(M);

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

void testChain(pimc::pimcConfigurations & configurations, int iChain, int expectedHead, int expectedTail)
{
    
    const auto chain = configurations.getChain(iChain);


    const auto & mask = configurations.getMask(); 

    ASSERT_EQ( chain.head , expectedHead   );
    ASSERT_EQ( chain.tail , expectedTail   );

    if (expectedHead > expectedTail)
    {
        for(int t=expectedTail+1;t<expectedHead;t++)
        {
            ASSERT_EQ(mask(t,iChain) , 1  );
        }
        for(int t=0;t<=expectedTail;t++)
        {
            ASSERT_EQ(mask(t,iChain) , 0  );
        }
        for(int t=expectedHead;t<configurations.nBeads();t++)
        {
            ASSERT_EQ(mask(t,iChain) , 0  );
        }

    }
    else 
    {
        for(int t=expectedTail;t<=expectedHead;t++)
        {
            ASSERT_EQ(mask(t,iChain) , 0  );
        }
        for(int t=0;t<expectedTail;t++)
        {
            ASSERT_EQ(mask(t,iChain) , 1  );
        }
        for(int t=expectedHead+1;t<configurations.nBeads();t++)
        {
            ASSERT_EQ(mask(t,iChain) , 1  );
        }

    }

}


TEST(configurations, worms)
{
    const int N = 10;
    const int M = 50;

    pimc::particleGroup groupA{ 0 , N-1, N + 1 , 1.0};

    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});
    
    auto & data = configurations.dataTensor();
    const auto & mask = configurations.getMask();

    data.setRandom();

    configurations.fillHeads();

    configurations.save("randomConfiguration");


    int time=10;
    int iChain=0;

    // test open
    configurations.setHead( iChain , M );

    ASSERT_EQ( configurations.heads().size() , 1 );
    ASSERT_EQ( configurations.tails().size() , 1 );
    
    ASSERT_EQ( configurations.heads()[0],iChain);
    ASSERT_EQ(configurations.tails()[0],iChain);


    ASSERT_TRUE( configurations.getChain(iChain).hasHead()     );
    ASSERT_TRUE(  configurations.getChain(iChain).hasTail() );

    // test close

    configurations.join(iChain,iChain);
    ASSERT_FALSE(configurations.getChain(iChain).hasHead() );
    ASSERT_FALSE(configurations.getChain(iChain).hasTail() );
    ASSERT_EQ( configurations.heads().size() , 0 );
    ASSERT_EQ( configurations.tails().size() , 0 );

    // test creation of a new head and join
     iChain=5;
    configurations.setHead(iChain,M);
    configurations.pushChain(0);
    configurations.setHead(N,time);
    configurations.setTail(N,-1);
    testChain(configurations, N, time, -1);
    testChain(configurations, iChain, M, -1); 
    configurations.join(iChain,N);

    ASSERT_EQ( configurations.heads().size() , 1 );
    ASSERT_EQ( configurations.tails().size() , 1 );


    ASSERT_EQ( configurations.tails()[0] , iChain );
    ASSERT_EQ( configurations.heads()[0] , N );
   


    // test creation of a new tail and join

    iChain=5;
    configurations.setTail(iChain,-1);
    int newChain=configurations.pushChain(0);
    configurations.setHead(newChain,M);
    configurations.setTail(newChain,time);
    testChain(configurations, newChain, M, time);
    testChain(configurations, iChain, M, -1); 
    configurations.join(newChain,iChain); 


    ASSERT_EQ( configurations.heads().size() , 1 );
    ASSERT_EQ( configurations.tails().size() , 1 );

    ASSERT_EQ( configurations.tails()[0] , newChain );
    ASSERT_EQ( configurations.heads()[0] , N );
    ASSERT_EQ(newChain,N+1);
    ASSERT_FALSE(configurations.getChain(N).hasTail() );
    ASSERT_FALSE(configurations.getChain(N+1).hasHead() );


    // test remove
    configurations.setTail(N,-1);

    ASSERT_EQ( configurations.tails().size() , 2 );

    configurations.setHead(N+1,0);
    
    ASSERT_EQ( configurations.tails()[1] , N );
    ASSERT_EQ( configurations.tails()[0] , N + 1 );
    
    configurations.removeChain(N); 
    configurations.removeChain(N);

    ASSERT_EQ( configurations.tails()[0] , iChain );
    ASSERT_EQ( configurations.heads()[0] , iChain );

}

TEST(run,free_harmonic_oscillator)
{   
    int N=10;
    int M=10;
    Real Beta = 1;
    

    pimc::geometryPBC_PIMC geo(300,300,300);

    Real timeStep = Beta/M;

    pimc::particleGroup groupA{ 0 , N-1, N - 1 , 1.0};
    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});

    configurations.dataTensor().setRandom();

    //configurations.join(0,1);
    //configurations.join(1,0);
    
    configurations.fillHeads();



    pimc::levyReconstructor reconstructor(M);

    pimc::levyMove freeMoves(reconstructor, 5);

    Real delta=0.1;

    pimc::translateMove translMove(delta,M*N);

    Real C = 1e-1;
    int l = 1;

    pimc::openMove openMove(C,l);
    pimc::closeMove closeMove(C,l);

    pimc::moveHead moveHeadMove(4);
    pimc::moveTail moveTailMove(l);

    pimc::swapMove swapMove(l,M);


    pimc::tableMoves table;
  
    table.push_back(& freeMoves,0.8,pimc::sector_t::offDiagonal,"levy");
    table.push_back(& freeMoves,0.8,pimc::sector_t::diagonal,"levy");

    table.push_back(& translMove,0.5,pimc::sector_t::diagonal,"translate");
    //table.push_back(& translMove,0.5,pimc::sector_t::offDiagonal,"translate");

    table.push_back(& openMove,0.2,pimc::sector_t::diagonal,"open");
    table.push_back(& closeMove,0.2,pimc::sector_t::offDiagonal,"close");

    table.push_back(& moveHeadMove,0.4,pimc::sector_t::offDiagonal,"moveHead");
    table.push_back(& moveTailMove,0.4,pimc::sector_t::offDiagonal,"moveTail");

    table.push_back(& swapMove,1.9,pimc::sector_t::offDiagonal,"swap");


    randomGenerator_t randG(19);

    pimc::kineticAction sT(timeStep, configurations.nChains() , M  , geo);

    #if DIMENSIONS == 3
    auto V = [](Real x, Real y , Real z) {return 0.5*(x*x + y*y + z*z) ;};
    #endif

    #if DIMENSIONS == 1
    auto V = pimc::makePotentialFunctor(
         [](Real x) {return 0.5*(x*x ) ;} ,
         [](Real x) {return 0.5*x ;} 
         );
    #endif

    pimc::potentialActionOneBody<decltype(V)> sV(timeStep,V ,geo);

    pimc::firstOrderAction S(&sT, & sV);
    int nTimes = 100000;
    int success = 0;
    int subSteps=10000;
    int correlationSteps=100;

    pimc::thermodynamicEnergyEstimator energyEstimator;
    Real e=0;
    Real e2=0;

    std::ofstream f;

    configurations.save("testRun/sample"+std::to_string(0));

    f.open("energy.dat");
    for (int i=0;i< nTimes ; i++)
    {
        Real eStep=0;
        int nMeasurements=0;

        for (int k=0;k< subSteps;k++)
        {
            
            for (int j=0;j<correlationSteps;j++)
            {
                bool accepted=table.attemptMove(configurations, S, randG);

                if (accepted)
                {success+=1;}
            }
            
            if (!configurations.isOpen() )
            {
                Real tmp=energyEstimator(configurations,S);
                nMeasurements++;
            
                eStep+=tmp;
            }
            else
            {
               
            }
            
        }

       
        f << i + 1 << "\t" << eStep/nMeasurements << std::endl ;

        //std::cout << e << std::endl;
        std::cout << "Energy: " << eStep/nMeasurements << std::endl;
        std::cout << "Acceptance ratio: " << success*1./((i+1)*subSteps*correlationSteps) << std::endl;

        table >> std::cout;

        configurations.save("testRun/sample"+std::to_string(i+1));
    }

    f.close();
    ASSERT_TRUE( (success*1./nTimes )> 0);
    std::cout << "END." << std::endl;
}