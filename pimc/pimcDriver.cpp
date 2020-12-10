#include "pimcDriver.h"
#include "geometryPMC.h"
#include "traits.h"
#include "tools.h"
#include "action.h"
#include "pimcConfigurations.h"
#include "moves.h"
#include "pimcObservables.h"
#include <filesystem>

namespace fs = std::filesystem;

namespace pimc
{

pimcDriver::pimcDriver(const json_t & j)
{
    std::vector<Real> lBox;
    lBox=j["particles"].get<std::vector<Real> >();

    Real beta = j["inverseTemperature"].get<Real>();
    nParticles = j["particles"].get<std::vector<int> >();

    nBeads= j["nBeads"].get<int>();

    timeStep = beta/nBeads;
    
    seed = j["seed"].get<int>();


    //int nChains = std::accumulate(nParticles.begin(),nParticles.end() , 0);
    int nChains = nParticles[0];
    


    #if DIMENSIONS == 1
    pimc::geometryPBC_PIMC geo(lBox[0]);
    #endif

    #if DIMENSIONS == 2
    pimc::geometryPBC_PIMC geo(lBox[0],lBox[1]);
    #endif
    
    #if DIMENSIONS == 3
    pimc::geometryPBC_PIMC geo(lBox[0],lBox[1],lBox[2]);
    #endif


    // build moves

    pimc::moveConstructor pimcMoveConstructor(nParticles,nBeads);
    pimcMoveConstructor.registerMove<pimc::levyMove>("levy");
    pimcMoveConstructor.registerMove<pimc::translateMove>("translate");
    pimcMoveConstructor.registerMove<pimc::openMove>("open");
    pimcMoveConstructor.registerMove<pimc::closeMove>("close");
    pimcMoveConstructor.registerMove<pimc::moveHead>("moveHead");
    pimcMoveConstructor.registerMove<pimc::moveTail>("moveTail");
    pimcMoveConstructor.registerMove<pimc::swapMove>("swap");
    
    tab = pimcMoveConstructor.createTable( j["movesTable"] );

    //auto swapMove = new pimc::swapMove(4,nParticles[0] );

    //tab.push_back(swapMove,1.9,pimc::sector_t::offDiagonal,"swap");

     stepsPerBlock = j["stepsPerBlock"].get<int>();
    nBlocks = j["nBlocks"].get<int>();
    correlationSteps = j["correlationSteps"].get<int>();


    




}

void pimcDriver::run()
{
    // build action 
    std::shared_ptr<action> sT= std::make_shared<kineticAction>(timeStep, nParticles[0] , nBeads  , geo);

     #if DIMENSIONS == 3
     auto V = pimc::makePotentialFunctor(
         [](Real x,Real y , Real z) {return 0.5*(x*x + y*y + z*z) ;} ,
         [](Real x,Real y, Real z) {return x  ;},
         [](Real x,Real y,Real z) {return y ;},
         [](Real x,Real y,Real z) {return z ;}
         );
    #endif

    #if DIMENSIONS == 1
    auto V = pimc::makePotentialFunctor(
         [](Real x) {return 0.5*(x*x ) ;} ,
         [](Real x) {return x ;} 
         );
    #endif


    std::shared_ptr<action> sV=std::make_shared<pimc::potentialActionOneBody<decltype(V)> >(timeStep,V ,geo);

    S = pimc::firstOrderAction(sT, sV);
    


   /*  pimc::levyMove freeMoves(5);

    Real delta=0.1;

    pimc::translateMove translMove(delta,(nBeads+1)*nParticles[0]);


   
     Real C = 1e-1;
    int l = 5;
    
    pimc::openMove openMove(C,l);
    pimc::closeMove closeMove(C,l);

    pimc::moveHead moveHeadMove(l);
    pimc::moveTail moveTailMove(l);

    pimc::swapMove swapMove(4,nParticles[0]);


     tab.push_back(& freeMoves,0.8,pimc::sector_t::offDiagonal,"levy");
    tab.push_back(& freeMoves,0.8,pimc::sector_t::diagonal,"levy");

    tab.push_back(& translMove,0.2,pimc::sector_t::diagonal,"translate");
    tab.push_back(& translMove,0.2,pimc::sector_t::offDiagonal,"translate");

    tab.push_back(& openMove,0.2,pimc::sector_t::diagonal,"open");
    tab.push_back(& closeMove,0.2,pimc::sector_t::offDiagonal,"close");

    tab.push_back(& moveHeadMove,0.4,pimc::sector_t::offDiagonal,"moveHead");
    tab.push_back(& moveTailMove,0.4,pimc::sector_t::offDiagonal,"moveTail");
    tab.push_back(& swapMove,1.9,pimc::sector_t::offDiagonal,"swap");
 */
    
     // build initial  configuration
    int N = nParticles[0];
    
    randomGenerator_t randG(seed);

    pimc::particleGroup groupA{ 0 , N-1, N - 1 , 1.0};

    pimc::pimcConfigurations configurations(nBeads, getDimensions() , {groupA});

    configurations.dataTensor().setRandom();
    configurations.fillHeads();


    pimc::thermodynamicEnergyEstimator energyEstimator;
    Real e=0;
    Real e2=0;

    std::ofstream f;

    if ( ! fs::exists("configurations") ) 
    { 
        fs::create_directory("configurations"); // create src folder
    }


    configurations.save("configurations/sample"+std::to_string(0));

    configurations.save("configurations/sample"+std::to_string(0));
    int success = 0;

    f.open("energy.dat");

    std::cout << "Start." << std::endl << std::flush;

    for (int i=0;i< nBlocks ; i++)
    {
        Real eStep=0;
        int nMeasurements=0;

        for (int k=0;k< stepsPerBlock;k++)
        {
            
            for (int j=0;j<correlationSteps;j++)
            {
                bool accepted=tab.attemptMove(configurations, S, randG);

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
        std::cout << "Acceptance ratio: " << success*1./((i+1)*stepsPerBlock*correlationSteps) << std::endl;

        tab >> std::cout;

        configurations.save("configurations/sample"+std::to_string(i+1));
    }
    
    f.close();
    std::cout << "END." << std::endl;

}


}