#include "gtest/gtest.h"
#include "geometryPMC.h"
#include "traits.h"
#include "tools.h"
#include "../pimc/action.h"
#include "../pimc/pimcConfigurations.h"


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

    //auto springDistances = geo.springDistances( timeConfigurations );

    //geo.springDistances(springDistances,timeConfigurations,0,2 , 0, 3  );

    
}

TEST(configurations, init)
{
    int N = 100;
    int M = 30;

    pimc::particleGroup groupA{ 0 , N-1, 1.0};

    pimc::pimcConfigurations configurations(M, N , getDimensions() , {groupA});

    auto & data = configurations.dataTensor();

    data.setRandom();

    pimc::geometryPBC_PIMC geo(10,10,10);


    Real timeStep = 1e-2;
    pimc::kineticAction sT(timeStep, N , M , geo);


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






}
