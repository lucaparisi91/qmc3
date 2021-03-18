#include "gtest/gtest.h"
#include "geometryPMC.h"
#include "traits.h"
#include "tools.h"
#include "../pimc/action.h"
#include "../pimc/pimcConfigurations.h"
#include "../pimc/moves.h"
#include "../pimc/pimcObservables.h"
#include "../pimc/hdf5IO.h"

#include <filesystem>

namespace fs = std::filesystem;


TEST(distances,updateSingleParticleDistances)
{
    pimc::geometryPBC_PIMC geo(10,10,10);
    int D = getDimensions();
    int N= 1;

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
    

    Eigen::Tensor<Real, 3> data(N,getDimensions(),T+1);
    data.setRandom();


    Eigen::Tensor<Real,3> springDifferences(T,getDimensions() , N );


    geo.updateSpringDifferences(  springDifferences, data, {0,T-1} , {0,N-1} );

    for (int t=0;t<T;t++)
    {
        for (int i=0;i<N;i++)
        {
            int d=0;
            ASSERT_NEAR( springDifferences(t,d,i) , data(i,d,t+1) - data(i,d,t) , 1e-4 );
        }
    }

    Eigen::Tensor<Real,3> potentialDifferences(N ,getDimensions() , T + 1 );

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

    pimc::particleGroup groupA{ 0 , N-1, N-1 , 1.0};
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

      #if DIMENSIONS == 3

      auto harmonicPotential = pimc::makeIsotropicPotentialFunctor(
         [](Real r) {return 0.5*(r*r) ;} ,
         [](Real r) {return r ;} 
         );
    #endif

      #if DIMENSIONS == 1
     auto harmonicPotential = pimc::makePotentialFunctor(
         [](Real x) {return 0.5*(x*x ) ;} ,
         [](Real x) {return x  ;},
         );
    #endif

    {
    Real x=0.1;
    Real y=0.23;
    Real z=0.56;
    
    ASSERT_NEAR(0.5*(x*x + y*y + z*z), harmonicPotential(x,y,z),1e-4);

    }

     pimc::potentialActionOneBody<decltype(harmonicPotential)> pot1(timeStep,harmonicPotential,geo);

     auto v = pot1.evaluate(configurations);
     Real vCheck=0;

     for (int t=1;t<=M-1;t++)
        for(int n=0;n<N;n++)
        {
            for(int d=0;d<getDimensions();d++)
            {
                vCheck+=0.5*data(n,d,t) * data(n,d,t) ;
            }
        }
    
      for(int n=0;n<N;n++)
        {
            for(int d=0;d<getDimensions();d++)
            {
                vCheck+=0.25*data(n,d,0) * data(n,d,0) ;
                vCheck+=0.25*data(n,d,M) * data(n,d,M) ;
            }
        }

        vCheck*=timeStep;
    ASSERT_NEAR(v,vCheck,1e-4);



}

#endif

#if DIMENSIONS == 3


TEST(moves,levy_reconstructor)
{
    int seed = 30;
    const int N = 50;
    const int M = 30;

    
    pimc::geometryPBC_PIMC geo(10,10,10);


    pimc::particleGroup groupA{ 0 , N-1, N ,  1.0};
    randomGenerator_t randG(seed+5);

    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});

    pimc::pimcConfigurations configurations2(M , getDimensions() , {groupA});

    Real timeStep=1e-3;

    std::array<int,2> timeSlice= {10,26};

    pimc::kineticAction sT(timeStep, configurations.nChains(), configurations.nBeads()  , geo);
     

    pimc::levyReconstructor levy(M);
    int iChain = 30;

    configurations.dataTensor().setRandom();
    configurations.fillHeads();

    configurations2=configurations;


    levy.apply(configurations,timeSlice,iChain,sT,randG);

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
        levy.apply(configurations,timeSlice,iChain,sT,randG);


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

     std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, N + 1 , M  , geo);
    
    #if DIMENSIONS == 1
    auto V = pimc::makePotentialFunctor(
         [](Real x) {return 0.5*(x*x ) ;} ,
         [](Real x) {return 0.5*x ;} 
         );    
    #endif

    #if DIMENSIONS == 3
    auto V = pimc::makeIsotropicPotentialFunctor(
         [](Real r) {return 0.5*(r*r) ;} ,
         [](Real r) {return r ;} 
         );    
    #endif


    std::shared_ptr<pimc::action> sV=std::make_shared<pimc::potentialActionOneBody<decltype(V)> >(timeStep,V ,geo);

    pimc::firstOrderAction S(sT,  sV);

    pimc::levyReconstructor levy(M);

    pimc::levyMove mover(20,0);
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
TEST(configurations, IO)
{
    const int N = 1000;
    const int M = 50;

    pimc::particleGroup groupA{ 0 , N-1, N , 1.0};

    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});

    auto & data = configurations.dataTensor();

    data.setRandom();

    configurations.fillHeads();

    configurations.join(100,120);
    configurations.join(120,100);

    configurations.setHead(10,M+1);

    std::string filename {"testConf.h5"} ;

    configurations.saveHDF5(filename);



     auto configurations2 = pimc::pimcConfigurations::loadHDF5(filename); 
    
    auto & data2 = configurations2.dataTensor();


    ASSERT_EQ(configurations.nBeads() , configurations2.nBeads() );
    ASSERT_EQ(configurations.nParticles() , configurations2.nParticles() );

    for(int t=0;t<M+1;t++)
        for (int i=0;i<N;i++)
            for(int d=0;d<getDimensions();d++)
                {
                   ASSERT_NEAR( data(i,d,t) , data2(i,d,t), 1e-5);
                } 

                
    for(int i=0;i<N;i++)
    {
        ASSERT_EQ( configurations2.getChain(i).prev , configurations.getChain(i).prev ); 
        ASSERT_EQ( configurations2.getChain(i).next , configurations.getChain(i).next ); 
        ASSERT_EQ( configurations2.getChain(i).head , configurations.getChain(i).head ); 
        ASSERT_EQ( configurations2.getChain(i).tail , configurations.getChain(i).tail );
        ASSERT_EQ( configurations2.getGroups()[0].tails[0] , configurations2.getGroups()[0].tails[0] );
    } 

    configurations.saveHDF5(filename);

}

TEST(observables, IO)
{
    const int N = 1000;
    const int M = 50;   
    Real Beta=1;
    Real timeStep=Beta/M;

    pimc::geometryPBC_PIMC geo(300,300,300);



    pimc::particleGroup groupA{ 0 , N-1, N , 1.0};

    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});

    auto & data = configurations.dataTensor();

    data.setRandom();

    configurations.fillHeads();

    
     std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);

    #if DIMENSIONS == 1
    auto V = pimc::makePotentialFunctor(
         [](Real x) {return 0.5*(x*x ) ;} ,
         [](Real x) {return 0.5*x ;} 
         );    
    #endif

    #if DIMENSIONS == 3

    auto V = pimc::makeIsotropicPotentialFunctor(
         [](Real r) {return 0.5*(r*r) ;} ,
         [](Real r) {return r ;} 
         );


    #endif


     std::shared_ptr<pimc::action> sV=std::make_shared<pimc::potentialActionOneBody<decltype(V)> >(timeStep,V ,geo);

    pimc::firstOrderAction S(sT,  sV);


    

    std::vector<std::shared_ptr<pimc::observable> > Os;
    {
        auto therm = std::make_shared<pimc::thermodynamicEnergyEstimator>();
        Os.push_back( std::make_shared<pimc::scalarObservable>(therm,"eT") );

    }


    for (auto & O : Os)
    {
        O->accumulate(configurations,S);
        O->out(0);
        O->clear();
        O->out(1);
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

    ASSERT_EQ( configurations.getGroup(0).heads.size() , 1 );
    ASSERT_EQ( configurations.getGroup(0).tails.size() , 1 );
    
    ASSERT_EQ( configurations.getGroup(0).heads[0],iChain);
    ASSERT_EQ(configurations.getGroup(0).tails[0],iChain);


    ASSERT_TRUE( configurations.getChain(iChain).hasHead()     );
    ASSERT_TRUE(  configurations.getChain(iChain).hasTail() );

    // test close

    configurations.join(iChain,iChain);
    ASSERT_FALSE(configurations.getChain(iChain).hasHead() );
    ASSERT_FALSE(configurations.getChain(iChain).hasTail() );
    ASSERT_EQ( configurations.getGroup(0).heads.size() , 0 );
    ASSERT_EQ( configurations.getGroup(0).tails.size() , 0 );

    // test creation of a new head and join
     iChain=5;
    configurations.setHead(iChain,M);
    configurations.pushChain(0);
    configurations.setHead(N,time);
    configurations.setTail(N,-1);
    testChain(configurations, N, time, -1);
    testChain(configurations, iChain, M, -1); 
    configurations.join(iChain,N);

    ASSERT_EQ( configurations.getGroup(0).heads.size() , 1 );
    ASSERT_EQ( configurations.getGroup(0).tails.size() , 1 );


    ASSERT_EQ( configurations.getGroup(0).tails[0] , iChain );
    ASSERT_EQ( configurations.getGroup(0).heads[0] , N );
    


    // test creation of a new tail and join

    iChain=5;
    configurations.setTail(iChain,-1);
    int newChain=configurations.pushChain(0);
    configurations.setHead(newChain,M);
    configurations.setTail(newChain,time);
    testChain(configurations, newChain, M, time);
    testChain(configurations, iChain, M, -1); 
    configurations.join(newChain,iChain); 


    ASSERT_EQ( configurations.getGroup(0).heads.size() , 1 );
    ASSERT_EQ( configurations.getGroup(0).tails.size() , 1 );

    ASSERT_EQ( configurations.getGroup(0).tails[0] , newChain );
    ASSERT_EQ( configurations.getGroup(0).heads[0] , N );
    ASSERT_EQ(newChain,N+1);
    ASSERT_FALSE(configurations.getChain(N).hasTail() );
    ASSERT_FALSE(configurations.getChain(N+1).hasHead() );


    // test remove
    configurations.setTail(N,-1);

    ASSERT_EQ( configurations.getGroup(0).tails.size() , 2 );

    configurations.setHead(N+1,0);
    
    ASSERT_EQ( configurations.getGroup(0).tails[1] , N );
    ASSERT_EQ( configurations.getGroup(0).tails[0] , N + 1 );
    
    configurations.removeChain(N); 
    configurations.removeChain(N);

    ASSERT_EQ( configurations.getGroup(0).tails[0] , iChain );
    ASSERT_EQ( configurations.getGroup(0).heads[0] , iChain );

}

TEST(action,twoBody)
{
    int N=100;
    int M=10;
    Real Beta = 1;

    Real timeStep=Beta/M;
    srand(11);

    pimc::geometryPBC_PIMC geo(300,300,300);

    pimc::particleGroup groupA{ 0 , N-1, N - 1 , 1.0};
    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});

    configurations.dataTensor().setRandom();
    configurations.fillHeads();


    auto & data = configurations.dataTensor();


    // Test on a rectangular interaction potential


    Real Rc=0.1;
    Real V0=2.;

    #if DIMENSIONS == 3


    auto V = pimc::makeIsotropicPotentialFunctor(
         [=](Real r) {return (r*r)<= Rc*Rc ? V0 :  0 ;} ,
         [](Real r) {return 0  ;}
         );
    #endif



    auto sV=pimc::potentialActionTwoBody<decltype(V)>(timeStep,N,M,V ,geo,0,0);

    int t0=0;
    int t1=M-1;
    int iChain=0;


    Real count=0;
    for(int t=t0;t<=t1;t++)
    {
         Real prefactor = 1;

        if ( (t==t0) or (t==t1) )
        {
            prefactor=0.5;
        }

            for(int j=0;j<N;j++)
            {
                
                Real dis=0;
                for(int d=0;d<getDimensions();d++)
                {
                    Real tmp=geo.difference(data(iChain,d,t)-data(j,d,t) ,d);
                    dis+=tmp*tmp;
                }

                if (dis<=Rc*Rc and j!= iChain)
                {
                    count+=prefactor;
                }
            }
    }

    auto currentAction = sV.evaluate(configurations,{t0,t1},iChain);
    ASSERT_NEAR(currentAction,count*V0*timeStep,1e-5);


    count=0;
    for(int t=0;t<=M-1;t++)
    {
        Real prefactor = 1;

        if ( (t==0) or (t==M-1) )
        {
            prefactor=0.5;
        }

        for(int i=0;i<N;i++)
            for(int j=0;j<i;j++)
            {
                Real dis=0;
                for(int d=0;d<getDimensions();d++)
                {
                    Real tmp=geo.difference(data(i,d,t)-data(j,d,t) ,d);
                    dis+=tmp*tmp;
                }

                if (dis<=Rc*Rc)
                {
                    
                    count+=prefactor;
                }
            }
    }

    auto totalAction = sV.evaluate(configurations);



    ASSERT_NEAR(totalAction,count*V0*timeStep,1e-5);

    Eigen::Tensor<Real,3> gradientBuffer(N,getDimensions(), M );
    Eigen::Tensor<Real,3> gradientBufferTest(N,getDimensions(), M );


    gradientBuffer.setConstant(0);
    gradientBufferTest.setConstant(0);
    

    sV.addGradient(configurations,{0,M-1},{0,N-1},gradientBuffer);

    Eigen::Tensor<Real,0> sumSquares = (gradientBuffer*gradientBuffer).sum();

    ASSERT_NEAR(sumSquares(0),0,1e-5);

    V0=1.;
    // test on a gaussian interaction potential
    Real alpha = 2.;


     auto V2 = pimc::makeIsotropicPotentialFunctor(
         [=](Real r) {return V0*exp(-alpha*(r*r));} ,
         [=](Real r) {return -2*r*V0*alpha*exp(-alpha*r*r)  ;}
          );
    

   
    auto sV2=pimc::potentialActionTwoBody<decltype(V2)>(timeStep,N,M,V2 ,geo,0,0);



    Real currentGaussV=0;
    for(int t=t0;t<=t1;t++)
    {
        Real prefactor = 1;
        if ( (t==0) or (t==(M-1)) )
        {
            prefactor=0.5;
        }

        for(int j=0;j<N;j++)
            {
                Real dis=0;
                for(int d=0;d<getDimensions();d++)
                {
                    Real tmp=geo.difference(data(iChain,d,t)-data(j,d,t) ,d);
                    dis+=tmp*tmp;
                }

                if (j != iChain)
                {
                    currentGaussV+=prefactor*V0*exp(-alpha*dis);
                }
                
            }
    }

    currentAction = sV2.evaluate(configurations,{t0,t1},iChain);
    ASSERT_NEAR(currentAction,currentGaussV*timeStep,1e-5);

    currentGaussV=0;
    for(int t=0;t<=M-1;t++)
    {
         Real prefactor = 1;
        if ( (t==0) or (t==M-1) )
        {
            prefactor=0.5;
        }
        for(int i=0;i<N;i++)
        for(int j=0;j<i;j++)
            {
                Real dis2=0;
                std::array<Real,3> diff;

                for(int d=0;d<getDimensions();d++)
                {
                    diff[d]=geo.difference(data(i,d,t)-data(j,d,t) ,d);
                    dis2+=diff[d]*diff[d];
                }

                currentGaussV+=prefactor*V0*exp(-alpha*dis2);

                for(int d=0;d<getDimensions();d++)
                {
                    gradientBufferTest(i,d,t)+=-2*alpha*diff[d]*V0*exp(-alpha*dis2)*timeStep;
                    gradientBufferTest(j,d,t)-=-2*alpha*diff[d]*V0*exp(-alpha*dis2)*timeStep;
                }

            }
    }

    totalAction = sV2.evaluate(configurations);
    ASSERT_NEAR(totalAction,currentGaussV*timeStep,1e-5);

    sV2.addGradient(configurations,{0,M-1},{0,N-1},gradientBuffer);

    Eigen::Tensor<Real,0> sumSquaresTest = (gradientBufferTest*gradientBufferTest).sum();

    sumSquares = (gradientBuffer*gradientBuffer).sum();

    ASSERT_NEAR(sumSquares(0),sumSquaresTest(0),1e-5);

    std::cout << sumSquares(0) << std::endl;
    std::cout << sumSquaresTest(0) << std::endl;
    
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

    pimc::levyMove freeMoves(5,0);

    Real delta=0.1;

    pimc::translateMove translMove(delta,(M+1)*N,0);


    Real C = 1e-1;
    int l = 3;

    
    pimc::openMove openMove(C,0,l);
    pimc::closeMove closeMove(C,0,l);

    pimc::moveHead moveHeadMove(l,0);
    pimc::moveTail moveTailMove(l,0);

    pimc::swapMove swapMove(l,N,0);

    pimc::tableMoves table;
    

    table.push_back(& freeMoves,0.8,pimc::sector_t::offDiagonal,"levy");
    table.push_back(& freeMoves,0.8,pimc::sector_t::diagonal,"levy");

    //table.push_back(& translMove,0.2,pimc::sector_t::diagonal,"translate");
    //table.push_back(& translMove,0.2,pimc::sector_t::offDiagonal,"translate");


    table.push_back(& openMove,0.2,pimc::sector_t::diagonal,"open");
    
    table.push_back(& closeMove,0.2,pimc::sector_t::offDiagonal,"close");

    table.push_back(& moveHeadMove,0.4,pimc::sector_t::offDiagonal,"moveHead");
    table.push_back(& moveTailMove,0.4,pimc::sector_t::offDiagonal,"moveTail");

    table.push_back(& swapMove,0.8,pimc::sector_t::offDiagonal,"swap");
    

    randomGenerator_t randG(368);

     std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);


     #if DIMENSIONS == 3
     auto V = pimc::makeIsotropicPotentialFunctor(
         [](Real r) {return 0.5*(r*r) ;} ,
         [](Real r) {return r  ;} );
    #endif

    #if DIMENSIONS == 1
    auto V = pimc::makePotentialFunctor(
         [](Real x) {return 0.5*(x*x ) ;} ,
         [](Real x) {return x ;} 
         );
    #endif

    Real R0=0.1;
    Real V0=1;



     #if DIMENSIONS == 3
  
      auto V2 = pimc::makeIsotropicPotentialFunctor(
         [=](Real r) {return V0*exp(-(r*r));} ,
         [=](Real r) {return -2*r*V0*exp(-r*r)  ;}
          );


    #endif

  

/* 
    #if DIMENSIONS == 1
    auto V2 = pimc::makePotentialFunctor(
         [=](Real x) {return  (x*x <=R0*R0) ? V0 : 0 ;} ,
         [](Real x) {return 0  ;}
         );        
    #endif */

    #if DIMENSIONS == 1
    auto V2 = pimc::makePotentialFunctor(
         [=](Real x) {return  exp(-x*x) ;} ,
         [](Real x) {return -2*x*exp(-x*x)  ;}
         );        
    #endif



    std::shared_ptr<pimc::action> sOneBody=std::make_shared<pimc::potentialActionOneBody<decltype(V)> >(timeStep,V ,geo);
    std::shared_ptr<pimc::action>  sV2B=std::make_shared<pimc::potentialActionTwoBody<decltype(V2)>  >(timeStep,N,M,V2 ,geo,0,0);    
    
    std::vector<std::shared_ptr<pimc::action> > Vs = {sOneBody,sV2B};
    
    
    std::shared_ptr<pimc::action>  sV = std::make_shared<pimc::sumAction>(Vs);

    pimc::firstOrderAction S(sT,  sV);
    
    int nTimes = 1000;
    int success = 0;
    int subSteps=1000;
    int correlationSteps=400;

   
    pimc::thermodynamicEnergyEstimator energyEstimator;

    pimc::virialEnergyEstimator viriralEnergy(N, M);

    Real e=0;
    Real e2=0;

    std::ofstream f;
    std::ofstream fV;


    if ( ! fs::exists("configurations") ) 
    { 
        fs::create_directory("configurations"); // create src folder
    }

    configurations.fillHeads();


    configurations.save("configurations/sample"+std::to_string(0),"pdb");

    f.open("energy.dat");
    fV.open("energyVirial.dat");

    for (int i=0;i< nTimes ; i++)
    {
        Real eStep=0,eVirialStep=0;
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
                Real tmp1=viriralEnergy(configurations,S);
            
                nMeasurements++;
            
                eStep+=tmp;
                eVirialStep+=tmp1;
            }
            else
            {
               
            }
            
        }

        f << i + 1 << "\t" << eStep/nMeasurements << std::endl ;
        fV << i + 1 << "\t" << eVirialStep/nMeasurements << std::endl ;

        //std::cout << e << std::endl;
        std::cout << "Energy: " << eStep/nMeasurements << std::endl;
        std::cout << "Acceptance ratio: " << success*1./((i+1)*subSteps*correlationSteps) << std::endl;

        table >> std::cout;

        configurations.save("configurations/sample"+std::to_string(i+1),"pdb");
    }

    f.close();
    ASSERT_TRUE( (success*1./nTimes )> 0);
    std::cout << "END." << std::endl;
}

TEST(run,free)
{   

    Real density=0.15884256651199277;
    int N=10;

    int M=10;
    Real Beta = 1;
    Real lBox = std::pow(N/density,1./3) ;
    

    pimc::geometryPBC_PIMC geo(lBox,lBox,lBox);


    Real timeStep = Beta/M;

    pimc::particleGroup groupA{ 0 , N-1, N - 1 , 1.0};
    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});

    auto & data = configurations.dataTensor(); 
    data.setRandom();

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        for(int d=0;d<getDimensions();d++)
        {
            data(i,d,j)=(data(i,d,j)-0.5 )*lBox;
        }
    }
    


    //configurations.join(0,1);
    //configurations.join(1,0);
    
    configurations.fillHeads();

    //configurations.dataTensor()(0,0,M)=configurations.dataTensor()(0,0,0) + 5*lBox;
    


    pimc::levyReconstructor reconstructor(M);

    pimc::levyMove freeMoves(M/3,0);


    Real delta=0.1;

    pimc::translateMove translMove(delta,(M+1)*N,0);

    Real C = 1e-4;
    int l = 1;


    pimc::openMove openMove(C,0,l);
    pimc::closeMove closeMove(C,0,l);

    pimc::moveHead moveHeadMove(l,0);
    pimc::moveTail moveTailMove(l,0);

    pimc::swapMove swapMove(4,N,0);

    pimc::tableMoves table;


    table.push_back(& freeMoves,0.8,pimc::sector_t::offDiagonal,"levy");
    table.push_back(& freeMoves,0.8,pimc::sector_t::diagonal,"levy");

    //table.push_back(& translMove,0.2,pimc::sector_t::diagonal,"translate");
    //table.push_back(& translMove,0.2,pimc::sector_t::offDiagonal,"translate");

    table.push_back(& openMove,0.1,pimc::sector_t::diagonal,"open");
    table.push_back(& closeMove,0.1,pimc::sector_t::offDiagonal,"close");

    table.push_back(& moveHeadMove,0.4,pimc::sector_t::offDiagonal,"moveHead");
    table.push_back(& moveTailMove,0.4,pimc::sector_t::offDiagonal,"moveTail");

    table.push_back(& swapMove,1.9,pimc::sector_t::offDiagonal,"swap");


    randomGenerator_t randG(368);

     std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);


     #if DIMENSIONS == 3
     auto V = pimc::makeIsotropicPotentialFunctor(
         [](Real r) {return 0 ;} ,
         [](Real r) {return 0  ;}
         );
    #endif


    #if DIMENSIONS == 1
    auto V = pimc::makePotentialFunctor(
         [](Real x) {return 0 ;} ,
         [](Real x) {return 0 ;} 
         );
    #endif

    Real alpha=1/std::pow (0.1*lBox,2);

    Real V0=100;

     #if DIMENSIONS == 3

      auto V2 = pimc::makeIsotropicPotentialFunctor(
         [=](Real r) {return V0*exp(-alpha*(r*r));} ,
         [=](Real r) {return -2*r*V0*alpha*exp(-alpha*r*r)  ;}
          );
    #endif


    #if DIMENSIONS == 1
    auto V2 = pimc::makePotentialFunctor(
         [=](Real x) {return  exp(-x*x) ;} ,
         [](Real x) {return -2*x*exp(-x*x)  ;}
         );        
    #endif



    std::shared_ptr<pimc::action> sOneBody=std::make_shared<pimc::potentialActionOneBody<decltype(V)> >(timeStep,V ,geo);
    std::shared_ptr<pimc::action>  sV2B=std::make_shared<pimc::potentialActionTwoBody<decltype(V2)>  >(timeStep,N,M,V2 ,geo,0,0);    
    
    std::vector<std::shared_ptr<pimc::action> > Vs = {sOneBody,sV2B};
    
    std::shared_ptr<pimc::action>  sV = std::make_shared<pimc::sumAction>(Vs);


    /*
    std::shared_ptr<pimc::action> sV=std::make_shared<pimc::potentialActionOneBody<decltype(V)> >(timeStep,V ,geo);
    */

    pimc::firstOrderAction S(sT,  sV);
    int nTimes = 1000;
    int success = 0;
    int subSteps=1000;
    int correlationSteps=100;

    pimc::thermodynamicEnergyEstimator energyEstimator;

    pimc::virialEnergyEstimator viriralEnergy(N, M);


    Real e=0;
    Real e2=0;

    std::ofstream f;
    std::ofstream fV;


    if ( ! fs::exists("configurations") ) 
    { 
        fs::create_directory("configurations"); // create src folder
    }


    configurations.save("configurations/sample"+std::to_string(0),"pdb");

    f.open("energy.dat");
    fV.open("energyVirial.dat");

    for (int i=0;i< nTimes ; i++)
    {
        Real eStep=0,eVirialStep=0;
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
                Real tmp1=viriralEnergy(configurations,S);
            
                nMeasurements++;
            
                eStep+=tmp;
                eVirialStep+=tmp1;
            }
            else
            {
               
            }
            
        }

        f << i + 1 << "\t" << eStep/nMeasurements << std::endl ;
        fV << i + 1 << "\t" << eVirialStep/nMeasurements << std::endl ;

        //std::cout << e << std::endl;
        std::cout << "Energy: " << eStep/nMeasurements << std::endl;
        std::cout << "Acceptance ratio: " << success*1./((i+1)*subSteps*correlationSteps) << std::endl;

        table >> std::cout;

        configurations.save("configurations/sample"+std::to_string(i+1),"pdb");
    }
    
    f.close();
    ASSERT_TRUE( (success*1./nTimes )> 0);
    std::cout << "END." << std::endl;
}

