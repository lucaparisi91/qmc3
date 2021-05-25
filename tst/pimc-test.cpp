#include "gtest/gtest.h"
#include "geometryPMC.h"
#include "traits.h"
#include "tools.h"
#include "../pimc/action.h"
#include "../pimc/pimcConfigurations.h"
#include "../pimc/moves.h"
#include "../pimc/pimcObservables.h"
#include "../pimc/hdf5IO.h"
#include "../pimc/toolsPimcTest.h"
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
        Os.push_back( std::make_shared<pimc::scalarObservable>(therm,std::string("eT") ) );

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
    const auto & tags = configurations.getTags();

    data.setRandom();

    configurations.fillHeads();

    configurations.save("randomConfiguration");


    int time=10;
    int iChain=0;

    // test open
    configurations.setHead( iChain , M );

    ASSERT_EQ( configurations.getGroupByChain(0).heads.size() , 1 );
    ASSERT_EQ( configurations.getGroupByChain(0).tails.size() , 1 );
    
    ASSERT_EQ( configurations.getGroupByChain(0).heads[0],iChain);
    ASSERT_EQ(configurations.getGroupByChain(0).tails[0],iChain);


    ASSERT_TRUE( configurations.getChain(iChain).hasHead()     );
    ASSERT_TRUE(  configurations.getChain(iChain).hasTail() );

    // test close

    configurations.join(iChain,iChain);
    ASSERT_FALSE(configurations.getChain(iChain).hasHead() );
    ASSERT_FALSE(configurations.getChain(iChain).hasTail() );
    ASSERT_EQ( configurations.getGroupByChain(0).heads.size() , 0 );
    ASSERT_EQ( configurations.getGroupByChain(0).tails.size() , 0 );

    // test creation of a new head and join
     iChain=5;
    configurations.setHead(iChain,M);
    configurations.pushChain(0);
    configurations.setHead(N,time);
    configurations.setTail(N,-1);
    testChain(configurations, N, time, -1);
    testChain(configurations, iChain, M, -1); 
    configurations.join(iChain,N);

    ASSERT_EQ( configurations.getGroupByChain(0).heads.size() , 1 );
    ASSERT_EQ( configurations.getGroupByChain(0).tails.size() , 1 );


    ASSERT_EQ( configurations.getGroupByChain(0).tails[0] , iChain );
    ASSERT_EQ( configurations.getGroupByChain(0).heads[0] , N );
    


    // test creation of a new tail and join

    iChain=5;
    configurations.setTail(iChain,-1);
    int newChain=configurations.pushChain(0);
    configurations.setHead(newChain,M);
    configurations.setTail(newChain,time);
    testChain(configurations, newChain, M, time);
    testChain(configurations, iChain, M, -1); 
    configurations.join(newChain,iChain); 


    ASSERT_EQ( configurations.getGroupByChain(0).heads.size() , 1 );
    ASSERT_EQ( configurations.getGroupByChain(0).tails.size() , 1 );

    ASSERT_EQ( configurations.getGroupByChain(0).tails[0] , newChain );
    ASSERT_EQ( configurations.getGroupByChain(0).heads[0] , N );
    ASSERT_EQ(newChain,N+1);
    ASSERT_FALSE(configurations.getChain(N).hasTail() );
    ASSERT_FALSE(configurations.getChain(N+1).hasHead() );


    // test remove
    configurations.setTail(N,-1);

    ASSERT_EQ( configurations.getGroupByChain(0).tails.size() , 2 );

    configurations.setHead(N+1,0);
    
    ASSERT_EQ( configurations.getGroupByChain(0).tails[1] , N );
    ASSERT_EQ( configurations.getGroupByChain(0).tails[0] , N + 1 );
    
    configurations.removeChain(N); 
    configurations.removeChain(N);

    ASSERT_EQ( configurations.getGroupByChain(0).tails[0] , iChain );
    ASSERT_EQ( configurations.getGroupByChain(0).heads[0] , iChain );

    pimc::particleGroup groupB{ 0 , N-1, N + 2 , 1.0};

    pimc::pimcConfigurations configurations2(M , getDimensions() , {groupB});
    
    configurations2.setEnsamble(pimc::ensamble_t::grandCanonical);
    int newHead=M/2;
    int newTail=M/2 ;
    int iOpen=N/2;

    configurations2.dataTensor().setRandom();
    configurations2.fillHeads();
    int iNext=configurations2.getChain(iOpen).next;

    configurations2.setHead(iOpen,newHead);
    auto iChainTail=configurations2.pushChain(0);
    configurations2.setTail(iChainTail,newTail);
    configurations2.join(iChainTail,iNext);

    ASSERT_EQ( configurations2.getGroupByChain(0).heads[0] , iOpen );
    ASSERT_EQ( configurations2.getChain(iChainTail).next , iNext );
    ASSERT_EQ( configurations2.getChain(iNext).prev , iChainTail );
    ASSERT_EQ( configurations2.getGroupByChain(0).tails[0] , iChainTail );
    ASSERT_EQ( configurations2.getGroupByChain(0).heads[0] , iOpen );


    testChain(configurations2,iOpen,newHead,-1);
    testChain(configurations2,iChainTail,M,newTail);

    configurations2.setHead(iChainTail,M);
    configurations2.removeChain(iOpen);

    testChain(configurations2,iOpen,M,newTail);


    ASSERT_EQ( configurations2.getGroupByChain(0).heads[0] , iOpen );
    ASSERT_EQ( configurations2.getGroupByChain(0).tails[0] , iOpen );


    configurations2.setTail(iOpen,-1);
    iChainTail=configurations2.pushChain(0);
    configurations2.setHead(iChainTail,M);
    configurations2.setTail(iChainTail,newTail);
    configurations2.join(iChainTail,iOpen);

    testChain(configurations2,iOpen,M,-1);
    testChain(configurations2,iChainTail,M,newTail);

}


TEST(action,oneBodyGrandCanonical)
{
    int N=100;
    int M=50;

    pimc::particleGroup group{ 0 , N-1, N + 2 , 1.0};

    pimc::pimcConfigurations configurations(M , getDimensions() , {group});

    configurations.setEnsamble(pimc::ensamble_t::grandCanonical);
    configurations.dataTensor().setRandom();

    int iOpen=N-1;
    int newHead=M/2;

    int iChainHead=configurations.pushChain(0);
    ASSERT_EQ(iChainHead,N);

    configurations.setHead(iOpen,M);
    configurations.setHead(iChainHead,newHead);
    configurations.setTail(iChainHead,-1);
    configurations.join(iOpen,iChainHead); 
    configurations.fillHeads();

    testChain(configurations,iChainHead,newHead,-1);
    
    pimc::geometryPBC_PIMC geo(300,300,300);

    Real timeStep=1e-1;

     auto V = pimc::makeIsotropicPotentialFunctor(
         [=](Real r) { return 0.5*r*r ;} ,
         [](Real r) {return r  ;}
         );


    const auto & mask = configurations.getTags();

    const auto & data = configurations.dataTensor();
    
    

    
     Real sum=0;
    for(int t=0;t<=M;t++)
    {
        Real prefactor = ((t ==0) or (t == M)) ? 0.5 : 1;
        for(int i=0;i<N;i++)
        {
            sum+=prefactor*V( TRUNCATE_D(data(i,0,t), data(i,1,t),data(i,2,t) ) );
        }
    }



     for(int t=0;t<=newHead;t++)
    {
        Real prefactor = (t ==0 or t == newHead) ? 0.5 : 1;
        int i=N;
        sum+=prefactor*V( TRUNCATE_D(data(i,0,t), data(i,1,t),data(i,2,t) ) );
    }


    auto sOneBody=std::make_shared<pimc::potentialActionOneBody<decltype(V)> >(timeStep,V ,geo);

    auto sumAction=sOneBody->evaluate(configurations);

    ASSERT_NEAR(sumAction,sum*timeStep,1e-5);
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

    auto V = pimc::makeIsotropicPotentialFunctor(
         [=](Real r) {return (r*r)<= Rc*Rc ? V0 :  0 ;} ,
         [](Real r) {return 0  ;}
         );


    auto sV=pimc::potentialActionTwoBody<decltype(V)>(timeStep,N,M,V ,geo,0,0);
    
    int t0=0;
    int t1=M-1;
    int iChain=0;

    Real count=0;
    for(int t=t0;t<=t1+1;t++)
    {
         Real prefactor = 1;

        if ( (t==t0) or (t==t1+1) )
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
    for(int t=0;t<=M;t++)
    {
        Real prefactor = 1;

        if ( (t==0) or (t==M) )
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



    currentAction = sV2.evaluate(configurations,{t0,t1-1},iChain);
    ASSERT_NEAR(currentAction,currentGaussV*timeStep,1e-5);

    currentGaussV=0;
    for(int t=0;t<=M;t++)
    {
         Real prefactor = 1;
        if ( (t==0) or (t==M) )
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

            }
    }

        for(int t=0;t<=M-1;t++)
    {
        
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
    int N=1;
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


    //table.push_back(& openMove,0.2,pimc::sector_t::diagonal,"open");
    
    table.push_back(& closeMove,0.2,pimc::sector_t::offDiagonal,"close");

    table.push_back(& moveHeadMove,0.4,pimc::sector_t::offDiagonal,"moveHead");
    table.push_back(& moveTailMove,0.4,pimc::sector_t::offDiagonal,"moveTail");

    table.push_back(& swapMove,0.8,pimc::sector_t::offDiagonal,"swap");
    

    randomGenerator_t randG(368);

     std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);


     auto V = pimc::makeIsotropicPotentialFunctor(
         [](Real r) {return 0.5*(r*r) ;} ,
         [](Real r) {return r  ;} );


    Real R0=0.1;
    Real V0=1;



  
      auto V2 = pimc::makeIsotropicPotentialFunctor(
         [=](Real r) {return V0*exp(-(r*r));} ,
         [=](Real r) {return -2*r*V0*exp(-r*r)  ;}
          );




    std::shared_ptr<pimc::action> sOneBody=std::make_shared<pimc::potentialActionOneBody<decltype(V)> >(timeStep,V ,geo);
    std::shared_ptr<pimc::action>  sV2B=std::make_shared<pimc::potentialActionTwoBody<decltype(V2)>  >(timeStep,N,M,V2 ,geo,0,0);    
    
    std::vector<std::shared_ptr<pimc::action> > Vs = {sOneBody};

    
    std::shared_ptr<pimc::action>  sV = std::make_shared<pimc::sumAction>(Vs);

    pimc::firstOrderAction S(sT,  sV);
    
    int nTimes = 1000;
    int success = 0;
    int subSteps=1000;
    int correlationSteps=10;

   
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


void accumulateBeadPosition(int i,std::array<Real,getDimensions()> & x, std::array<Real,getDimensions()> & x2, const pimc::configurations_t & configurations, const pimc::firstOrderAction & S)
    {
        const auto & data = configurations.dataTensor();
        for(int d=0;d<getDimensions();d++)
            {
                    x[d]+=data( 0 ,d, i) ;
                    x2[d]+=std::pow(data(0,d,i),2);

                }
    }


Real accumulateX2( const pimc::configurations_t & configurations, const pimc::firstOrderAction & S)
{
    Real x2=0;
    const auto & data = configurations.dataTensor();
    const auto & group = configurations.getGroups()[0];    

    for (int t=0;t<configurations.nBeads();t++)
        for (int i=group.iStart;i<=group.iEnd;i++)
            {
                for(int d=0;d<getDimensions();d++)
                {
                    x2+=data( i ,d, t)*data(i,d,t);
                }

            }
    

    return x2/configurations.nBeads();
            
}

Real accumulateX2SingleDirection(int iChain, const pimc::configurations_t & configurations, int direction, bool &  isCyclic, int t0 , int t1)
    {
        const auto & data = configurations.dataTensor();
        Real l2=0;

        int iCurrentChain=iChain;

        if (iCurrentChain < 0)
        {
            return 0;
        }
        do 
        {
            const auto & chain = configurations.getChain(iCurrentChain);

            for (int t=std::max(chain.tail + 1,t0);t<=std::min(chain.head , t1 );t++)
            {
                Real prefactor = (t == chain.tail+1) or (t == chain.head) ? 0.5 : 1;
                for(int d=0;d<getDimensions();d++)
                {
                        l2+=prefactor*data( iCurrentChain ,d, t)*data( iCurrentChain ,d, t);
                }
            }

            if (direction == 1)
            {
                iCurrentChain=chain.next;
            }
            else if (direction == -1)
            {
                iCurrentChain=chain.prev;
            }

        }
        while ( (iCurrentChain!= -1) and (iCurrentChain != iChain) );


        isCyclic= iCurrentChain == - 1 ? false : true;
        
        return l2;

    }

Real accumulateAverageLengthSquareSingleDirection(int iChain, const pimc::configurations_t & configurations, int direction, bool &  isCyclic, int t0 , int t1)
    {
        const auto & data = configurations.dataTensor();
        Real l2=0;

        int iCurrentChain=iChain;

        if (iCurrentChain < 0)
        {
            return 0;
        }
        do 
        {
            const auto & chain = configurations.getChain(iCurrentChain);

            for (int t=std::max(chain.tail + 1,t0);t<=std::min(chain.head - 1, t1 );t++)
            {
                for(int d=0;d<getDimensions();d++)
                {
                        l2+=std::pow(data( iCurrentChain ,d, t+1) - data(iCurrentChain,d,t),2);
                }
            }

            if (direction == 1)
            {
                iCurrentChain=chain.next;
            }
            else if (direction == -1)
            {
                iCurrentChain=chain.prev;
            }

        }
        while ( (iCurrentChain!= -1) and (iCurrentChain != iChain) );


        isCyclic= iCurrentChain == - 1 ? false : true;
        
        return l2;

    }


Real accumulateAverageLengthSquare(int iChain, const pimc::configurations_t & configurations, int t0, int t1)
{
    bool isCyclic;

    Real l2=accumulateAverageLengthSquareSingleDirection(iChain, configurations,+1,isCyclic,t0,t1);

    if (not isCyclic)
    {
        
        l2+=accumulateAverageLengthSquareSingleDirection(configurations.getChain(iChain).prev, configurations,-1,isCyclic,t0,t1);

        assert(isCyclic == false);

    }

    return l2;
}


Real accumulateX2(int iChain, const pimc::configurations_t & configurations, int t0, int t1)
{
    bool isCyclic;

    Real l2=accumulateX2SingleDirection(iChain, configurations,+1,isCyclic,t0,t1);

    if (not isCyclic)
    {
        
        l2+=accumulateX2SingleDirection(configurations.getChain(iChain).prev, configurations,-1,isCyclic,t0,t1);

        assert(isCyclic == false);

    }

    return l2;
}




Real accumulateAverageLengthSquare(int iChain, const pimc::configurations_t & configurations)
{
    int t0=0;
    int t1=configurations.nBeads()-1;



    return accumulateAverageLengthSquare( iChain, configurations,0,t1);
}


class configurationsTest : public ::testing::Test {
protected:
    configurationsTest() {
        open=NULL;
        close=NULL;
    }

     void SetRandom()
    {
        std::uniform_real_distribution<double> uniformDistribution(0.0,1.0);

        auto & data=configurations.dataTensor();
        for (int t=0;t<data.dimensions()[2];t++)
        for (int i=0;i<data.dimensions()[0];i++)
            for  (int d=0;d<getDimensions();d++)
            {
                data(i,d,t)=uniformDistribution(randG);
            }
        configurations.fillHeads();
    
    }

    void SetUp( int N_ , int M_ , Real Beta_) {
        N=N_;
        M=M_;
        Beta=Beta_;

        int seed= 356;
        int buffer=2;

        int nChains=N + buffer;

        timeStep=Beta/M;
        pimc::particleGroup groupA{ 0 , N-1, nChains -1 , 1.0};


        configurations=pimc::pimcConfigurations(M , getDimensions() , {groupA});

        geo = pimc::geometryPBC_PIMC(300,300,300);

        randG=randomGenerator_t(seed);

        SetRandom();

    }

    void SetSeed(int seed)
    {
        randG=randomGenerator_t(seed);
    }

    void SetUpFreeParticleAction()
    {
         std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);

    std::shared_ptr<pimc::action> sV= std::make_shared<pimc::nullPotentialAction>(timeStep  , geo);


    
    S= pimc::firstOrderAction(sT,  sV);



    
    }

    void SetUpNonInteractingHarmonicAction()
    {
         std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);




     auto V = pimc::makeIsotropicPotentialFunctor(
         [](Real r) {return 0.5*(r*r) ;} ,
         [](Real r) {return r  ;} );



    std::shared_ptr<pimc::action> sOneBody=std::make_shared<pimc::potentialActionOneBody<decltype(V)> >(timeStep,V ,geo);
    
    
    std::vector<std::shared_ptr<pimc::action> > Vs = {sOneBody};

    std::shared_ptr<pimc::action>  sV = std::make_shared<pimc::sumAction>(Vs);

    S= pimc::firstOrderAction(sT,  sV);

    }
    

    void SetGrandCanonicalEnsamble(Real chemicalPotential)
    {
        configurations.setEnsamble(pimc::ensamble_t::grandCanonical);
        configurations.setChemicalPotential(chemicalPotential);
    }

   

    void addOpenCloseFixedSegment(Real C , int l , int t0)
    {
        open=new pimc::openMove(C, 0, l );
        close=new pimc::closeMove(C, 0, l );

        open->setStartingBead(t0);
        close->setStartingBead(t0);

        open->setLengthCut(l);
        close->setLengthCut(l);

        tab.push_back(open,1,pimc::sector_t::diagonal);
        tab.push_back(close, 1, pimc::sector_t::offDiagonal);
    }



    
    template<class T>
    void accumulate(int nBurns,int nTrials,const T & f, int correlationSteps=1, pimc::sector_t sector= pimc::sector_t::diagonal)
    {
        for(int k=0;k<nBurns;k++)
        {
            bool accept= tab.attemptMove(configurations,S,randG);    
        }

        for(int k=0;k<nTrials;k++)
        {
            for (int kk=0;kk<correlationSteps;kk++)
            {
                 bool accept= tab.attemptMove(configurations,S,randG);

            }


            nMoves+=1;
            if ( configurations.isOpen() )
            {
                    

                    nOpen+=1;

                    if (sector == pimc::sector_t::offDiagonal)
                    {
                        f(configurations,S);
                    }
                    
                
            }
            else
            {
                nClosed+=1;

                  if (sector == pimc::sector_t::diagonal)
                    {
                        f(configurations,S);
                    }


            }

           

        }


    }



    virtual void TearDown() {
        if (open != NULL )
        {
            delete open;
        }

        if (close != NULL )
        {
            delete close;
        }

    }

    void resetCounters()
    {
        nMoves=0;
        nOpen=0;
        nClosed=0;
    }


    int N;
    int M;
    Real Beta;
    Real timeStep;
    pimc::pimcConfigurations configurations;
    pimc::firstOrderAction S;
    randomGenerator_t randG;
    pimc::geometryPBC_PIMC geo;
    pimc::tableMoves tab;

    pimc::openMove *open;
    pimc::closeMove *close;

    Real nOpen=0;
    Real nClosed=0;
    Real nMoves=0;

    
};



auto meanBeadFixedLengths(int iChainBegin , int iChainEnd,  int t0, int t1, int i, Real timeStep, const pimc::configurations_t & configurations)
{
    const auto & data = configurations.dataTensor();



    std::array<Real,getDimensions()> meanExpected;
    
    Real D = 0.5;
    Real mass = 1;
    int M= configurations.nBeads();

    std::array<Real,3> difference;
    for (int d=0;d<getDimensions();d++)
    {
        int l1 = (i - t0) > 0 ? i - t0 : i - t0 + M;
        int l2 = (t1 - i ) > 0 ? t1 - i  : t1 - i  + M;


        meanExpected[d]=(data(iChainBegin,d,t0)/ l1 + 
        data(iChainEnd,d,t1)/l2 )/(1./l1 + 1./l2 );

    }

    return meanExpected;
}

auto varianceBeadFixedLengths(int iChain , int t0, int t1, int i, Real timeStep, const pimc::configurations_t & configurations)
{

    std::array<Real,getDimensions()> varianceExpected;
    const auto & data = configurations.dataTensor();

    Real D = 0.5;
    Real mass = 1;

    int M= configurations.nBeads();

    std::array<Real,3> difference;
    for (int d=0;d<getDimensions();d++)
    {
        difference[d]= data(iChain,d,t0) - data(iChain,d,t1);

        int l1 = (i - t0) > 0 ? i - t0 : i - t0 + M;
        int l2 = (t1 - i ) > 0 ? t1 - i  : t1 - i  + M;


       varianceExpected[d]=1./(1./l1 + 1./l2 )* 2 * D * timeStep / mass;




    }

    return varianceExpected;
}


TEST_F(configurationsTest,openCloseGrandCanonical_distributionReconstructedChain)
{
    Real C=1;
    int t0=80;
    int l =40;

    SetUp(1,100,1);
    SetGrandCanonicalEnsamble(0);
    SetUpFreeParticleAction();
    addOpenCloseFixedSegment(C , l , t0);

    int t1=(t0 + l)%M;
    int i = (t0 + 20)%M;

    Real nTrials = 1000000;
    Real nBurns = 100;

    std::array<Real,getDimensions() > x{0,0,0};
    std::array<Real ,getDimensions() > x2{0,0,0};

    auto & data = configurations.dataTensor();
    

    pimc::configurations_t configurationsInitial=configurations;

    const auto & dataInitial = configurations.dataTensor();


    accumulate(nBurns,nTrials, [&](const pimc::configurations_t & confs, const pimc::firstOrderAction & S){accumulateBeadPosition(i,x,x2,confs,S) ;}  );


    Real openRatio = nOpen/nMoves;
    
    auto  meanExpected = meanBeadFixedLengths(0,0,t0,t1,i,timeStep,configurations);

    auto  varianceExpected = varianceBeadFixedLengths(0,t0,t1,i,timeStep,configurations);
    

  
    Real D = 0.5;
    Real mass = 1;


    std::array<Real,3> difference;
    for (int d=0;d<getDimensions();d++)
    {
          difference[d]= data(0,d,t0) - data(0,d,t1);
        x[d]/=nClosed;
        x2[d]/=nClosed;

    }

    Real expectedOpenRatio= 1/( 1 + exp(pimc::freeParticleLogProbability(difference,S.getTimeStep()*l,mass))/C );

    std::cout << "Open ratio error : " << std::abs(openRatio - expectedOpenRatio)/expectedOpenRatio <<  std::endl;

    EXPECT_NEAR((openRatio - expectedOpenRatio)/expectedOpenRatio , 0 , 1e-2);


    std::cout << "Mean error: " << std::abs((x[0]  -  meanExpected[0])/meanExpected[0])<< std::endl;
    EXPECT_NEAR( (x[0]  -  meanExpected[0])/meanExpected[0] , 0 , 1e-2);

    std::cout << "var error : " << std::abs(x2[0] - x[0]*x[0] - varianceExpected[0] )/varianceExpected[0] << std::endl;

    EXPECT_NEAR(( x2[0] - x[0]*x[0] - varianceExpected[0] )/varianceExpected[0] , 0 , 3e-2);

}


TEST_F(configurationsTest,openChain)
{
    Real C=1;
    int nBeads=10;

    SetUp(2,nBeads,1);
    //SetGrandCanonicalEnsamble(0);
    SetUpFreeParticleAction();
    
    SetRandom();


    SetSeed(time(NULL));


    
    int t0=10;
    int l = nBeads/3;

    pimc::levyMove levy(l,0);
    pimc::moveHead moveHeadMove(l,0);
    pimc::moveTail moveTailMove(l,0);



    configurations.setTail(0,-1);
    configurations.setHead(1,M);
    configurations.join(0,1);

    configurations.fillHeads();

    tab.push_back(&levy,0.9,pimc::sector_t::offDiagonal);
    tab.push_back(&moveHeadMove,0.9,pimc::sector_t::offDiagonal);
    tab.push_back(&moveTailMove,0.9,pimc::sector_t::offDiagonal);


    Real l2=0;
    Real l2Var=0;
    Real l2Error=0;

    int nTrials = 1000000;

    accumulate(1000,nTrials, [&](const pimc::configurations_t & confs, const pimc::firstOrderAction & S){Real tmp=accumulateAverageLengthSquare(0,confs) ;l2+=tmp;l2Var+=tmp*tmp;}  ,10,pimc::sector_t::offDiagonal);

    l2/=nOpen;
    l2Var/=nOpen;

    l2Error = std::sqrt((l2Var - l2*l2)/nOpen);

    EXPECT_NEAR(l2 , 3* 2, 2*l2Error);


}


TEST_F(configurationsTest,closedChain_free)
{
    Real C=1;
    int nBeads=10;

    SetUp(2,nBeads,1);
    //SetGrandCanonicalEnsamble(0);
    SetUpFreeParticleAction();
    
    SetRandom();


    SetSeed(time(NULL));

    
    int t0=10;
    int l = nBeads/3;

    pimc::levyMove levy(l,0);




    configurations.setTail(0,-1);
    configurations.setHead(1,M);
    configurations.join(0,1);
    configurations.join(1,0);


    //configurations.fillHeads();

    tab.push_back(&levy,0.9,pimc::sector_t::diagonal);

    Real l2=0;
    Real l2Var=0;
    Real l2Error=0;

    int nTrials = 100000;

    accumulate(1000,nTrials, [&](const pimc::configurations_t & confs, const pimc::firstOrderAction & S){Real tmp=accumulateAverageLengthSquare(0,confs) ;l2+=tmp;l2Var+=tmp*tmp;}  ,100,pimc::sector_t::diagonal);

    l2/=nClosed;
    l2Var/=nClosed;

    l2Error = std::sqrt((l2Var - l2*l2)/nClosed);

    
    EXPECT_NEAR(l2 , 3* (2*M - 1) * timeStep , 2*l2Error);


}


TEST_F(configurationsTest,closedChain_harmonic)
{
    Real C=1;
    int nBeads=10;

    SetUp(1,nBeads,1);
    //SetGrandCanonicalEnsamble(0);
    SetUpNonInteractingHarmonicAction();

    SetRandom();


    SetSeed(time(NULL));

    
    int t0=10;
    int l = nBeads/3;

    pimc::levyMove levy(l,0);

/* 

    configurations.setTail(0,-1);
    configurations.setHead(1,M);
    configurations.join(0,1);
    configurations.join(1,0);


    //configurations.fillHeads();
 */
    tab.push_back(&levy,0.9,pimc::sector_t::diagonal);

    Real l2=0;
    Real l2Var=0;
    Real l2Error=0;
    int i=4;

    int nTrials = 100000000;
    std::array<Real,3> x2 {0,0,0};
    std::array<Real,3> x {0,0,0};


    accumulate(1000,nTrials, [&](const pimc::configurations_t & confs, const pimc::firstOrderAction & S){Real tmp=accumulateAverageLengthSquare(0,confs) ;l2+=tmp;l2Var+=tmp*tmp;
    
    accumulateBeadPosition(i, x, x2, confs, S);
    
    }  ,100,pimc::sector_t::diagonal);

    for (int d=0;d<getDimensions();d++)
    {
        x[d]/=nClosed;
        x2[d]/=nClosed;
    }


    l2/=nClosed;
    l2Var/=nClosed;

    l2Error = std::sqrt((l2Var - l2*l2)/nClosed);

    //EXPECT_NEAR(l2 , 3* (2*M - 1) * timeStep , 2*l2Error);

    std::cout << l2 << " " << l2Error << std::endl;

    for (int d=0;d<getDimensions();d++)
    {
        std::cout << x[d] << " " << x2[d] << std::endl;

    }

}


TEST_F(configurationsTest,openClosedChain_harmonic)
{
    Real C=1;
    int nBeads=10;

    SetUp(1,nBeads,1);
    SetGrandCanonicalEnsamble(0);
    SetUpNonInteractingHarmonicAction();

    SetSeed(356);

    SetRandom();
    
    int t0=9;
    int l = 4;

    pimc::levyMove levy(l,0);

    pimc::openMove open(C, 0, l );
    pimc::closeMove close(C, 0, l );

    pimc::moveHead moveHeadMove(l,0);
    pimc::moveTail moveTailMove(l,0);
    

    open.setStartingBead(t0);
    close.setStartingBead(t0);

    open.setLengthCut(l);
    close.setLengthCut(l);

    tab.push_back(&levy,0.9,pimc::sector_t::diagonal);
    tab.push_back(&levy,0.9,pimc::sector_t::offDiagonal);
    tab.push_back(&moveHeadMove,0.1,pimc::sector_t::offDiagonal);
    tab.push_back(&moveTailMove,0.1,pimc::sector_t::offDiagonal);
    

    tab.push_back(&close,0.9,pimc::sector_t::offDiagonal);
    tab.push_back(&open,0.9,pimc::sector_t::diagonal);
    


    //configurations.setHeadTail(0,M,-1);

    //configurations.setTail(0,-1);
    //configurations.setHead(1,M);
    //configurations.join(0,1);
    //configurations.join(1,0);


    //configurations.fillHeads();


    //Real l2=0;
    //Real l2Var=0;
    //Real l2Error=0;
   // int iNewHead=4;

    int nTrials = 10000;
    
    std::ofstream l2Out,x2Out,openFractionOut;

    l2Out.open("l2.dat");
    x2Out.open("x2.dat");
    openFractionOut.open("openFraction.dat");
    

    for (int iBlock=0;iBlock < 1000000 ; iBlock++)
    {
        Real l2=0;
        Real l2Var=0;
        Real x2=0;

        accumulate(0,nTrials, [&](
            
            const pimc::configurations_t & confs, const pimc::firstOrderAction & S){

                auto tmp=accumulateAverageLengthSquare(0,confs) ;l2+=tmp;
                x2+=accumulateX2(0, configurations, 0, M);

                const auto & tags = confs.getTags();



                
            
            } 
                ,20,pimc::sector_t::diagonal
        );


        int nMeasures= nClosed;

        l2Out << iBlock << "\t" << l2/nMeasures << std::endl;
        x2Out << iBlock << "\t" << x2/nMeasures << std::endl;
        openFractionOut << iBlock <<  "\t" << nOpen/(nOpen + nClosed) << std::endl;


        resetCounters();



    }


    l2Out.close();
    x2Out.close();
    openFractionOut.close();


/* 

    l2/=nClosed;
    l2Var/=nClosed;

    l2Error = std::sqrt((l2Var - l2*l2)/nClosed);

    //EXPECT_NEAR(l2 , 3* (2*M - 1) * timeStep , 2*l2Error);

    std::cout << l2 << " " << l2Error << std::endl;

    for (int d=0;d<getDimensions();d++)
    {
        std::cout << x[d] << " " << x2[d] << std::endl;

    } */



}


TEST_F(configurationsTest,openClosedChain_free)
{
    Real C=1e-1;
    int nBeads=10;

    SetUp(1,nBeads,1);
    SetGrandCanonicalEnsamble(0);
    SetUpFreeParticleAction();

    SetSeed(time(NULL));

    SetRandom();
    
    int t0=3;
    int l = nBeads/3;

    pimc::levyMove levy(l,0);

    pimc::openMove open(C, 0, l );
    pimc::closeMove close(C, 0, l );

    pimc::moveHead moveHeadMove(l,0);
    pimc::moveTail moveTailMove(l,0);

    open.setStartingBead(t0);
    close.setStartingBead(t0);

    open.setLengthCut(l);
    close.setLengthCut(l);


    tab.push_back(&levy,0.9,pimc::sector_t::diagonal);
    tab.push_back(&levy,0.9,pimc::sector_t::offDiagonal);
    tab.push_back(&moveHeadMove,0.1,pimc::sector_t::offDiagonal);
    tab.push_back(&moveTailMove,0.1,pimc::sector_t::offDiagonal);


    tab.push_back(&close,0.9,pimc::sector_t::offDiagonal);
    tab.push_back(&open,0.9,pimc::sector_t::diagonal);
    


    //configurations.setHeadTail(0,M,-1);

    //configurations.setTail(0,-1);
    //configurations.setHead(1,M);
    //configurations.join(0,1);
    //configurations.join(1,0);


    //configurations.fillHeads();


    //Real l2=0;
    //Real l2Var=0;
    //Real l2Error=0;
    int iNewHead=4;

    int nTrials = 10000;
    
    std::ofstream l2Out,x2Out;

    l2Out.open("l2.dat");
    x2Out.open("x2.dat");

    for (int iBlock=0;iBlock < 1000000 ; iBlock++)
    {
        Real l2=0;
        Real l2Var=0;
        Real x2=0;

        accumulate(0,nTrials, [&](
            
            const pimc::configurations_t & confs, const pimc::firstOrderAction & S){

                auto tmp=accumulateAverageLengthSquare(0,confs) ;l2+=tmp;
                //x2+=accumulateX2(0, configurations, 0, M);

            } 
                ,20,pimc::sector_t::offDiagonal
        );


        int nMeasures= nOpen;

        std::cout << nOpen/(nClosed + nOpen) << std::endl;

        l2Out << iBlock << "\t" << l2/nMeasures << std::endl;
        //x2Out << iBlock << "\t" << x2/nMeasures << std::endl;
        
        resetCounters();

    }
    
    l2Out.close();
    x2Out.close();


/* 

    l2/=nClosed;
    l2Var/=nClosed;

    l2Error = std::sqrt((l2Var - l2*l2)/nClosed);

    //EXPECT_NEAR(l2 , 3* (2*M - 1) * timeStep , 2*l2Error);

    std::cout << l2 << " " << l2Error << std::endl;

    for (int d=0;d<getDimensions();d++)
    {
        std::cout << x[d] << " " << x2[d] << std::endl;

    } */

}


TEST_F(configurationsTest,openChain_harmonic)
{
    Real C=1;
    int nBeads=10;

    SetUp(2,nBeads,1);
    SetGrandCanonicalEnsamble(0);
    SetUpNonInteractingHarmonicAction();

    SetRandom();


    SetSeed(time(NULL));

    
    int t0=10;
    int l = 3;
    int iHead=6;
    int iTail=5;

    pimc::levyMove levy(l,0);
    pimc::moveHead moveHeadMove(l,0);
    pimc::moveTail moveTailMove(l,0);



    configurations.setHeadTail(1,M,iTail);
    configurations.setHeadTail(0,iHead,-1);
    configurations.join(1,0);
    configurations.fillHeads();

    



    tab.push_back(&levy,0.9,pimc::sector_t::offDiagonal);
    tab.push_back(&moveHeadMove,0.2,pimc::sector_t::offDiagonal);
    tab.push_back(&moveTailMove,0.2,pimc::sector_t::offDiagonal);



    int nTrials = 1000;

    std::ofstream x2Out,l2Out;


    l2Out.open("l2.dat");
    x2Out.open("x2.dat");

    for (int iBlock=0;iBlock < 1000000 ; iBlock++)
    {
            Real l2=0;
            Real x2=0;

            accumulate(
                1000,nTrials, [&](
                
                const pimc::configurations_t & confs, const pimc::firstOrderAction & S){Real tmp=accumulateAverageLengthSquare(0,confs) ;l2+=tmp;
                x2+=accumulateX2(0, configurations, 0, M);

                }
                ,100,pimc::sector_t::offDiagonal
                    );
            
            l2Out << iBlock << " " << l2/nOpen << std::endl;
            x2Out << iBlock << " " << x2/nOpen << std::endl;
        
            resetCounters();

    
    }
    
    l2Out.open("l2Out.dat");
    x2Out.open("x2.dat");



/*     l2/=nOpen;
    l2Var/=nOpen;

    l2Error = std::sqrt((l2Var - l2*l2)/nOpen);

    //EXPECT_NEAR(l2 , 3* (2*M - 1) * timeStep , 2*l2Error);

    std::cout << l2 << " " << l2Error << std::endl;

    for (int d=0;d<getDimensions();d++)
    {
        std::cout << x[d] << " " << x2[d] << std::endl;

    } */

}






TEST_F(configurationsTest,openClosedChain)
{
    Real C=1;
    int nBeads=10;

    SetUp(1,nBeads,1);
    SetGrandCanonicalEnsamble(0);
    SetUpFreeParticleAction();
    
    SetRandom();


    SetSeed(time(NULL));

    
    int t0=1;
    int l = 3;


    pimc::levyMove levy(l,0);

    pimc::openMove open(C, 0, l );
    pimc::closeMove close(C, 0, l );

    pimc::moveHead moveHeadMove(l,0);
    pimc::moveTail moveTailMove(l,0);



    tab.push_back(&levy,0.9,pimc::sector_t::diagonal);
    tab.push_back(&levy,0.9,pimc::sector_t::offDiagonal);
    tab.push_back(&moveHeadMove,0.1,pimc::sector_t::offDiagonal);
    tab.push_back(&moveTailMove,0.1,pimc::sector_t::offDiagonal);


    tab.push_back(&close,2,pimc::sector_t::offDiagonal);
    tab.push_back(&open,2,pimc::sector_t::diagonal);

    open.setStartingBead(t0);
    close.setStartingBead(t0);

    open.setLengthCut(l);
    close.setLengthCut(l);


    Real l2=0;
    Real l2Var=0;
    Real l2Error=0;

    int nTrials = 1000000;

    accumulate(1000,nTrials, [&](const pimc::configurations_t & confs, const pimc::firstOrderAction & S){Real tmp=accumulateAverageLengthSquare(0,confs) ;l2+=tmp;l2Var+=tmp*tmp;}  ,100,pimc::sector_t::diagonal);

    l2/=nClosed;
    l2Var/=nClosed;

    l2Error = std::sqrt((l2Var - l2*l2)/nClosed);

    EXPECT_NEAR(l2 , 3* (M - 1) * timeStep , 3*l2Error);


    resetCounters();


    l2=0;
    l2Var=0;

    accumulate(1000,nTrials, [&](const pimc::configurations_t & confs, const pimc::firstOrderAction & S){Real tmp=accumulateAverageLengthSquare(0,confs) ;l2+=tmp;l2Var+=tmp*tmp;}  ,100,pimc::sector_t::offDiagonal);

    l2/=nOpen;
    l2Var/=nOpen;

    l2Error = std::sqrt((l2Var - l2*l2)/nOpen);


    EXPECT_NEAR(l2 , 3* ( M - l ) * timeStep , 3*l2Error);

    
}




TEST_F(configurationsTest,advanceRecedeGrandCanonical_distributionReconstructedChain)
{
    Real C=1;
    int l =40;
    int nBeads= 100;

    SetUp(1,nBeads,1);
    SetGrandCanonicalEnsamble(0);
    SetUpFreeParticleAction();
    
    int iHead=60;

    configurations.setHeadTail(0,iHead,-1);

    pimc::advanceHead advanceMove(l,0);
    pimc::recedeHead recedeMove(l,0);

    advanceMove.setFixedLength();
    recedeMove.setFixedLength();

    int nTrails = 1000000;

    std::array<Real ,getDimensions()> x{0,0,0}; 
    int nMeasurements = 0;
    int i = iHead + l <= M ? iHead + l  : iHead + l - M;
    int iChainNewHead = iHead + l <= M ? 0  : 1;

    const auto & data = configurations.dataTensor();

    SetSeed(56);

    for (int n=0;n<nTrails;n++)
    {
        int currentHead = configurations.getChain(0).head;
        
        if (currentHead == iHead)
        {
            advanceMove.attemptMove(configurations,S,randG);
        }
        else
        {
            for (int d=0;d<getDimensions();d++)
            {
                x[d]+=data(iChainNewHead,d,i);
            }
            nMeasurements+=1;
            recedeMove.attemptMove(configurations,S,randG);
        }

    }

    ASSERT_NEAR( nMeasurements * 1./nTrails, 0.5 , 1e-3 ) ;

    for (int d=0;d<getDimensions();d++)
    {
        x[d]/=nMeasurements;
        ASSERT_NEAR ( data(0,d,iHead)  , x[d] , 1e-2 );

    }



}

TEST_F(configurationsTest,swapGrandCanonical_distributionReconstructedChain)
{
    Real C=1;
    int l =40;
    int nBeads= 100;

    SetUp(2,nBeads,1);
    SetGrandCanonicalEnsamble(0);
    SetUpFreeParticleAction();
    

    int iHead=30;

    configurations.setHeadTail(0,iHead,-1);

    pimc::swapMove swap(l,2,0);

    swap.setFixedLength();

    SetSeed(58);

    int nTrails = 1000000;
    int nBurns = 1000;
    int nMeasurements=0;

    const auto & data = configurations.dataTensor();

    pimc::configurations_t configurationsInitial(configurations);
    const auto & dataInitial = configurations.dataTensor();

    

    std::array<Real ,getDimensions()> x{0,0,0}; 
    int i = iHead + 10;

    auto meanExpected=meanBeadFixedLengths(0 ,  1 , iHead, iHead + l, i , timeStep,configurations);


    for (int n=0;n<nBurns;n++)
    {
        bool accept=swap.attemptMove(configurations,S,randG);
    }
    for (int n=0;n<nTrails;n++)
    {
        bool accept=swap.attemptMove(configurations,S,randG);

        assert( configurations.getChain(0).hasHead() or configurations.getChain(1).hasHead( ));


        if ( not configurations.getChain(0).hasHead()  )
        {
            
            for (int d=0;d<getDimensions();d++)
            {
                x[d]+=data(0,d,i);
            }
            nMeasurements+=1;
        }

    }

    for (int d=0;d<getDimensions();d++)
            {
                x[d]/=nMeasurements;
            }
    
    

    for (int d=0;d<getDimensions();d++)
    {

       

        int iChainHead = configurations.getGroups()[0].heads[0];
        int iChainCont = iChainHead == 0 ? 1 : 0;

            for (int t=0;t<=iHead;t++)
            {
                assert(data(0,d,t) == dataInitial(0,d,t)    );
                assert(data(1,d,t) == dataInitial(1,d,t)    );
            }

            for (int t=iHead+l;t<=M;t++)
            {
                assert(data(iChainCont,d,t) == dataInitial(1,d,t)    );
            }
         ASSERT_NEAR(x[d],meanExpected[d],1e-2);

           }
        
        
    
}




TEST(moves,openCloseGrandCanonical)
{   
    int N=1;
    int M=100;
    Real Beta = 1;

    int seed=812; // 356
    int buffer=2;

    int nChains=N + buffer;

    Real timeStep=Beta/M;
    pimc::particleGroup groupA{ 0 , N-1, nChains -1 , 1.0};

    std::default_random_engine generator;
    std::uniform_real_distribution<double> uniformDistribution(0.0,1.0);

    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});


    configurations.setEnsamble(pimc::ensamble_t::grandCanonical);
    configurations.setChemicalPotential(0);

    pimc::geometryPBC_PIMC geo(300,300,300);

    
    std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);

    std::shared_ptr<pimc::action> sV= std::make_shared<pimc::nullPotentialAction>(timeStep  , geo);


    randomGenerator_t randG(seed);

   

    pimc::firstOrderAction S(sT,  sV);

    auto & data=configurations.dataTensor();

    for (int t=0;t<data.dimensions()[2];t++)
        for (int i=0;i<data.dimensions()[0];i++)
            for  (int d=0;d<getDimensions();d++)
            {
                data(i,d,t)=uniformDistribution(randG);
            }
    
    Real C = 1;
    int l=30;
    int t0=6;
    int t1=t0 + l;
    int i = t0 + 20;



    pimc::tableMoves tab;

    pimc::openMove open(C, 0, l );
    pimc::closeMove close(C, 0, l );

    open.setStartingBead(t0);
    close.setStartingBead(t0);

    open.setLengthCut(l);
    close.setLengthCut(l);

    tab.push_back(& open,1,pimc::sector_t::diagonal);
    tab.push_back(& close, 1, pimc::sector_t::offDiagonal);


    configurations.fillHeads();

    Real nTrials = 100000;
    Real nBurns = 100;

    Real nOpen=0;

    std::array<Real,getDimensions() > x{0,0,0};
    std::array<Real ,getDimensions() > x2{0,0,0};


    for(int k=0;k<nBurns;k++)
    {
        bool accept= tab.attemptMove(configurations,S,randG);    
    }
    
    for(int k=0;k<nTrials;k++)
    {

        if ( configurations.isOpen() )
        {
            nOpen+=1;
        }
        else
        {
                for(int d=0;d<getDimensions();d++)
                {
                    x[d]+=data( 0 ,d, i) ;
                    x2[d]+=std::pow(data(0,d,i),2);

                }
            
            
        }

        bool accept= tab.attemptMove(configurations,S,randG);

    }

    Real nClosed = nTrials - nOpen;

    Real openRatio = nOpen/nTrials;
    
    std::array<Real,getDimensions()> meanExpected;
    std::array<Real,getDimensions()> varianceExpected;

    Real D = 0.5;
    Real mass = 1;

    std::array<Real,3> difference;
    for (int d=0;d<getDimensions();d++)
    {
        difference[d]= data(0,d,t0) - data(0,d,t1);

        varianceExpected[d]=1./(1./(i - t0) + 1./(t1 - i) )* 2 * D * timeStep / mass;

        meanExpected[d]=(data(0,d,t0)/ (i-t0) + 
        data(0,d,t1)/( t1 - i) )/(1./(i - t0) + 1./(t1 - i) );

        x[d]/=nClosed;
        x2[d]/=nClosed;


    }


    Real expectedOpenRatio= 1/( 1 + exp(pimc::freeParticleLogProbability(difference,S.getTimeStep()*l,mass))/C );

    std::cout << "Open ratio error : " << std::abs(openRatio - expectedOpenRatio)/expectedOpenRatio <<  std::endl;

    ASSERT_NEAR((openRatio - expectedOpenRatio)/expectedOpenRatio , 0 , 1e-2);


    std::cout << "Mean error: " << std::abs((x[0]  -  meanExpected[0])/meanExpected[0])<< std::endl;
    ASSERT_NEAR( (x[0]  -  meanExpected[0])/meanExpected[0] , 0 , 1e-2);

    std::cout << "var error : " << std::abs(x2[0] - x[0]*x[0] - varianceExpected[0] )/varianceExpected[0] << std::endl;

    ASSERT_NEAR(( x2[0] - x[0]*x[0] - varianceExpected[0] )/varianceExpected[0] , 0 , 1e-2);








}

TEST(run,free_harmonic_oscillator_grandCanonical)
{   
    int N=1;
    int M=10;
    Real Beta = 1;

    int seed= time(NULL);
    int buffer=200;

    int nChains=N + buffer;

    std::srand((unsigned int) seed);
    randomGenerator_t randG(seed);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> uniformDistribution(0.0,1.0);

    pimc::geometryPBC_PIMC geo(300,300,300);

    Real timeStep = Beta/M;

    pimc::particleGroup groupA{ 0 , N-1, nChains -1 , 1.0};

    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});

    //configurations.setEnsamble(pimc::ensamble_t::grandCanonical);
    //configurations.setChemicalPotential(3/2. * 0);


    auto & data=configurations.dataTensor();

    for (int t=0;t<data.dimensions()[2];t++)
        for (int i=0;i<data.dimensions()[0];i++)
            for  (int d=0;d<getDimensions();d++)
            {
                data(i,d,t)=uniformDistribution(randG);
            }


    std::cout << configurations.dataTensor().sum();
    //configurations.join(0,1);
    //configurations.join(1,0);
    
    configurations.fillHeads();

    pimc::levyReconstructor reconstructor(M);

    pimc::levyMove freeMoves(5,0);

    Real delta=0.1;

    //pimc::translateMove translMove(delta,(M+1)*N,0);

    Real C = 1e-3;
    int l = 4;
    
    pimc::openMove openMove(C,0,l);
    pimc::closeMove closeMove(C,0,l);

    pimc::moveHead moveHeadMove(l,0);
    pimc::moveTail moveTailMove(l,0);


    pimc::advanceHead advanceHeadMove(l,0);
    pimc::recedeHead recedeHeadMove(l,0);

    pimc::createWorm addMove(C,0,l,1);
    pimc::deleteWorm removeMove(C,0,l,1);

    pimc::swapMove swapMove(l,N,0);

    pimc::tableMoves table;

    table.push_back(& freeMoves,0.8,pimc::sector_t::offDiagonal,"levy");
    table.push_back(& freeMoves,0.8,pimc::sector_t::diagonal,"levy");

    //table.push_back(& translMove,0.2,pimc::sector_t::diagonal,"translate");
    //table.push_back(& translMove,0.2,pimc::sector_t::offDiagonal,"translate");

    table.push_back(& openMove,0.2,pimc::sector_t::diagonal,"open");
    table.push_back(& closeMove,0.2,pimc::sector_t::offDiagonal,"close");


    //table.push_back(& addMove,0.2,pimc::sector_t::diagonal,"addMove");
    //table.push_back(& removeMove,0.2,pimc::sector_t::offDiagonal,"removeMove");


    //table.push_back(& moveHeadMove,0.4,pimc::sector_t::offDiagonal,"moveHead");
    //table.push_back(& moveTailMove,0.4,pimc::sector_t::offDiagonal,"moveTail");

    //table.push_back(& advanceHeadMove,0.9,pimc::sector_t::offDiagonal,"advanceHead");

    //table.push_back(& recedeHeadMove,0.9,pimc::sector_t::offDiagonal,"recedeHead");

    //table.push_back(& swapMove,0.8,pimc::sector_t::offDiagonal,"swap");


     std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);


     auto V = pimc::makeIsotropicPotentialFunctor(
         [](Real r) {return 0.5*(r*r) ;} ,
         [](Real r) {return r  ;} );


    Real R0=0.1;
    Real V0=1;



  
      auto V2 = pimc::makeIsotropicPotentialFunctor(
         [=](Real r) {return V0*exp(-(r*r));} ,
         [=](Real r) {return -2*r*V0*exp(-r*r)  ;}
          );




    std::shared_ptr<pimc::action> sOneBody=std::make_shared<pimc::potentialActionOneBody<decltype(V)> >(timeStep,V ,geo);
    std::shared_ptr<pimc::action>  sV2B=std::make_shared<pimc::potentialActionTwoBody<decltype(V2)>  >(timeStep,nChains,M,V2 ,geo,0,0);    
    
    std::vector<std::shared_ptr<pimc::action> > Vs = {sOneBody};
    
    
    std::shared_ptr<pimc::action>  sV = std::make_shared<pimc::sumAction>(Vs);

    pimc::firstOrderAction S(sT,  sV);
    
    int nTimes = 100000;
    int success = 0;
    int subSteps=10000;
    int correlationSteps=10;
   
    pimc::thermodynamicEnergyEstimator energyEstimator;

    pimc::virialEnergyEstimator viriralEnergy(nChains, M);

    pimc::particleNumber particleEstimator(0);

    Real e=0;
    Real e2=0;
    Real sumN=0;


    std::ofstream f;
    std::ofstream fV;
    std::ofstream fN;
    

    if ( ! fs::exists("configurations") ) 
    { 
        fs::create_directory("configurations"); // create src folder
    }

    configurations.fillHeads();

    configurations.save("configurations/sample"+std::to_string(0),"pdb");

    f.open("energy.dat",std::ios_base::app);
    fV.open("energyVirial.dat",std::ios_base::app);
    fN.open("N.dat",std::ios_base::app);

    int iSuspicious=0;
    int iOpen=0;


    for (int i=0;i< nTimes ; i++)
    {
        Real eStep=0,eVirialStep=0,nStep=0;
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

                Real tmp2=particleEstimator(configurations,S);
            
                nMeasurements++;

                /* if (tmp < -6)
                {
                    std::string name="configurations/sampleSusp"+std::to_string(iSuspicious);
                    configurations.save(name,"pdb"); 
                    iSuspicious++;             
                } */
                eStep+=tmp;
                eVirialStep+=tmp1;
                nStep+=tmp2;
            }
            else
            {
                //configurations.save("configurations/sampleOpen"+std::to_string(iOpen),"pdb");
                //iOpen++;
            }
            
        }

        f << i + 1 << "\t" << eStep/nMeasurements << std::endl ;
        fV << i + 1 << "\t" << eVirialStep/nMeasurements << std::endl ;

        fN << i + 1 << "\t" << nStep/nMeasurements << std::endl ;

        //std::cout << e << std::endl;
        std::cout << "Energy: " << eStep/nMeasurements << std::endl;
        std::cout << "Acceptance ratio: " << success*1./((i+1)*subSteps*correlationSteps) << std::endl;
        std::cout << "N: " << configurations.nParticles() << std::endl;

        table >> std::cout;

        //configurations.save("configurations/sample"+std::to_string(i+1),"pdb");
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

    Real C = 1e-1;
    int l = 4;


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


   
     auto V = pimc::makeIsotropicPotentialFunctor(
         [](Real r) {return 0 ;} ,
         [](Real r) {return 0  ;}
         );

    Real alpha=1/std::pow (0.1*lBox,2);

    Real V0=100;


      auto V2 = pimc::makeIsotropicPotentialFunctor(
         [=](Real r) {return V0*exp(-alpha*(r*r));} ,
         [=](Real r) {return -2*r*V0*alpha*exp(-alpha*r*r)  ;}
          );


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

