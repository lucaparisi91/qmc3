#include "gtest/gtest.h"
#include "geometryPMC.h"
#include "traits.h"
#include "tools.h"
#include "../pimc/action.h"
#include "../pimc/pimcConfigurations.h"
#include "../pimc/moves.h"

TEST(distances,updateSingleParticleDistances)
{
    pimc::geometryPBC_PIMC geo(10,10,10);
    int D = getDimensions();
    int N= 5;

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
    const int N = 100;
    const int M = 30;

    pimc::particleGroup groupA{ 0 , N-1, 1.0};

    pimc::pimcConfigurations configurations(M, N , getDimensions() , {groupA});

    auto & data = configurations.dataTensor();

    data.setRandom();

    pimc::geometryPBC_PIMC geo(10,10,10);


    Real timeStep = 1e-2;
    pimc::kineticAction sT(timeStep, N , M , {groupA} , geo);

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

     pimc::potentialActionOneBody<decltype(harmonicPotential)> pot1(timeStep,harmonicPotential, {groupA},geo);

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

    pimc::particleGroup groupA{ 0 , N-1, 1.0};

    pimc::pimcConfigurations configurations(M, N , getDimensions() , {groupA});

    auto & data = configurations.dataTensor();

    data.setRandom();

    pimc::geometryPBC_PIMC geo(10,10,10);


    Real timeStep = 1e-2;
    pimc::kineticAction sT(timeStep, N , M , {groupA} , geo);

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
    const int N = 100;
    const int M = 30;

    pimc::particleGroup groupA{ 0 , N-1, 1.0};
    pimc::pimcConfigurations configurations(M, N , getDimensions() , {groupA});

    pimc::pimcConfigurations configurations2(M, N , getDimensions() , {groupA});

    Real timeStep=1e-3;

    std::array<int,2> timeSlice= {10,26};


    pimc::levyReconstructor levy(seed,timeStep);
    int iChain = 30;
    configurations2=configurations;

    levy.apply(configurations,configurations,iChain,timeSlice);

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
            ASSERT_FALSE( abs( data2(i,d,t) - data(i,d,t) ) < 1e-4);
        }

        };



    }


    levy.apply(configurations,configurations, 8,{25,40} );
    configurations2=configurations;

}

TEST(moves,levy)
{
    int seed = 30;
    const int N = 100;
    const int M = 30;

    pimc::particleGroup groupA{ 0 , N-1, 1.0};
    pimc::pimcConfigurations configurations(M, N , getDimensions() , {groupA});
    pimc::geometryPBC_PIMC geo(10,10,10);


    pimc::pimcConfigurations configurations2(M, N , getDimensions() , {groupA});
    Real timeStep=1e-3;
    std::array<int,2> timeSlice= {10,26};

    pimc::kineticAction sT(timeStep, N , M , {groupA} , geo);


    auto V = [](Real x, Real y , Real z) {return x*x + y*y + z*z ;};


    pimc::potentialActionOneBody<decltype(V)> sV(timeStep,V , {groupA},geo);

    pimc::firstOrderAction S(&sT, & sV);


    pimc::levyReconstructor levy(seed,timeStep);
    pimc::levyMove mover(levy,20);

    mover.attemptMove(configurations,S);
    







    

}