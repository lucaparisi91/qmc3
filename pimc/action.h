#include "unsupported/Eigen/CXX11/Tensor"
#include "geometryPMC.h"
#include "pimcConfigurations.h"
#include "tools.h"

namespace pimc
{


    template<class functor_t>
    Real reduceOnSpringDistances(const functor_t & V,const Eigen::Tensor<Real,3> & tn, std::array<int ,2 > timeRange, std::array<int, 2>  particleRange)
    {
        Real sum=0;

         for (size_t i=particleRange[0];i<=particleRange[1];i++ )
            for(int t=timeRange[0];t<=timeRange[1] ; t++ )
                {
                    sum+=V( tn( t,0, i  ) , tn(t,1,i) , tn(t,2,i)  ) ;
                }

        return sum;
    };

    template<class functor_t>
    Real reduceOnPositions(const functor_t & V,const Eigen::Tensor<Real,3> & tn, std::array<int ,2 > timeRange, std::array<int, 2>  particleRange)
    {
        Real sum=0;

         for (size_t i=particleRange[0];i<=particleRange[1];i++ )
            for(int t=timeRange[0];t<=timeRange[1] ; t++ )
                {
                    sum+=V( tn( i,0, t  ) , tn(i,1,t) , tn(i,2,t)  ) ;
                }

        return sum;
    };









class harmonicPotential
{
    inline Real operator() ( Real  x , Real y , Real z )
    {
        return 0.5*(x*x + y*y + z*z);
    }

};


using pimcConfigurations_t = pimcConfigurations;


class kineticAction
{
    public:
    kineticAction(Real tau_, int nChains_ , int nBeads_, const geometryPBC_PIMC & geo_) : geo(geo_),tau(tau_),nChains(nChains_),nBeads(nBeads_),distancesBuffer(nBeads_, getDimensions( ) , nChains_),
    D(0.5)  {}

    Real evaluate( pimcConfigurations_t & configurations , std::array<int,2> timeSlices , std::array<int,2> chainRange ); // evaluates the kinetic action to recalculate from beads(iParticle, timeSliceStart: timeSliceEnd, uses in internal buffer for computing derivatives)
    

    Real evaluate(pimcConfigurations_t & configurations); // evaluates the full action


    private:
    geometryPBC_PIMC geo;
    int nChains;
    int nBeads;


    Eigen::Tensor<Real, 3> distancesBuffer;
    Real tau;
    Real D;

};


template<class functor_t>
class potentialActionOneBody
{
    public:

    potentialActionOneBody(Real tau_, functor_t V_, geometryPBC_PIMC geo_) : tau(tau_),V(V_),geo(geo_) {}
    
    Real evaluate(pimcConfigurations_t & configurations, std::array<int,2> timeRange, std::array<int,2>  particleRange  ) 
    {
        auto & data = configurations.dataTensor();

        auto sum = reduceOnPositions(V, data, timeRange, particleRange);

        return sum;
    }



    private:
    functor_t V;
    Real tau;
    geometryPBC_PIMC geo;
    

};





template<class functor_t>
class potentialActionTwoBody
{
    public:


    potentialActionTwoBody(Real tau_, int nChains_  , int nBeads_, functor_t V_, geometryPBC_PIMC geo_) : tau(tau_),V(V_),geo(geo_),nChains(nChains_),nBeads(nBeads_) ,
    bufferDistances(nChains  , getDimensions() , nBeads  ) {}

    
    Real evaluate(pimcConfigurations_t & configurations, std::array<int,2> timeRange, int iChain, std::array<int,2>  particleRange  ) 
    {
        auto & data = configurations.dataTensor();

        geo.updateEqualTimeDifferences(bufferDistances, data, timeRange , iChain, particleRange) ;

        auto sum=reduceOnDifferences(bufferDistances,timeRange,particleRange);


        return sum;
    }

    Real evaluate(const pimcConfigurations_t & configurations)
    {
        Real sum=0;
        for(int iChain=0 ; iChain < nChains; iChain++ )
        {
            sum+=evaluate(configurations, {0,nBeads-1} , iChain, {0, nChains-1});
        }
        return sum/2.;
    };
    




    private:


    Real reduceOnDifferences(const Eigen::Tensor<Real,3> & tn, std::array<int,2> timeRange, std::array<int,2> particleRange  ) const 
    {
         Real sum=0;
        for(int t=timeRange[0];t<=timeRange[1] ; t++)
            for (size_t i=particleRange[0];i<=particleRange[1];i++ )
                {
                    sum+=V( tn(i,0,t ) , tn(i,1,t) , tn(i,2,t)  ) ;
                }

        return sum;

    }
    functor_t V;
    Real tau;
    Eigen::Tensor<Real,3> bufferDistances;
    geometryPBC_PIMC geo;
    int nChains,nBeads;

};


/* class action
{
    public:
    action(real_t tau);
    real_t operator()(configurations_t & conf);

    private:
    potential_t oneBodyPotential;
    real_t tau;

}; */

};

