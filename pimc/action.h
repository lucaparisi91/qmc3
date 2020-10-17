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
    Real reduceOnSpringDistances(const functor_t & V,const Eigen::Tensor<Real,3> & tn, std::array<int ,2 > timeRange, std::array<int, 2>  particleRange, const mask & mask)
    {
        Real sum=0;

         for (size_t i=particleRange[0];i<=particleRange[1];i++ )
            for(int t=timeRange[0];t<=timeRange[1] ; t++ )
                {
                    sum+= ( mask(i,t)==0 ? 0 : V( tn( t,0, i  ) , tn(t,1,i) , tn(t,2,i)  ) );
                }
        return sum;
    };

    template<class functor_t>
    Real reduceOnPositions(const functor_t & V,const Eigen::Tensor<Real,3> & tn, std::array<int ,2 > timeRange, std::array<int, 2>  particleRange)
    {
        Real sum=0;

        
            for(int t=timeRange[0];t<=timeRange[1] ; t++ )
             for (size_t i=particleRange[0];i<=particleRange[1];i++ )
                {
                    sum+=V( tn( i,0, t  ) , tn(i,1,t) , tn(i,2,t)  ) ;
                }

        return sum;
    };


    template<class functor_t>
    Real reduceOnPositions(const functor_t & V,const Eigen::Tensor<Real,3> & tn, std::array<int ,2 > timeRange, std::array<int, 2>  particleRange, const mask & mask)
    {
        Real sum=0;

         
            for(int t=timeRange[0];t<=timeRange[1] ; t++ )
              for (size_t i=particleRange[0];i<=particleRange[1];i++ )
                {
                    sum+= (mask(i,t) == 0 ? 0 :  V( tn( i,0, t  ) , tn(i,1,t) , tn(i,2,t)  ) ) ;
                }

        return sum;
    };

using pimcConfigurations_t = pimcConfigurations;


class action
{
    public:
    virtual Real evaluate(configurations_t & configurations, std::array<int,2> timeRange, int iParticle)=0;

    //virtual Real evaluate(const pimcConfigurations_t & configurations, std::array<int,2> timeRange, int iChain1 , int iChain2  )=0;

    virtual Real evaluate( configurations_t & pimcConfigurations)=0; // evaluates the whole action

    virtual Real evaluate( configurations_t & configurations , std::array< std::array<int,2> ,2> timeRanges, int iParticle)
    {
        return evaluate(configurations,timeRanges[0] , iParticle ) + evaluate(configurations,timeRanges[1],iParticle);
        
    }

};







class kineticAction : public action
{
    public:
    using action::evaluate;


    kineticAction(Real tau_, int nChains_ , int nBeads_, std::vector<particleGroup>  groups, const geometryPBC_PIMC & geo_
    ) : geo(geo_),tau(tau_),nChains(nChains_),nBeads(nBeads_),distancesBuffer(nBeads_, getDimensions( ) , nChains_),
    D(0.5) , _particleGroups(groups) {}


    Real evaluate( pimcConfigurations_t & configurations , std::array<int,2> timeSlices , std::array<int,2> chainRange ); // evaluates the kinetic action to recalculate from beads(iParticle, timeSliceStart: timeSliceEnd, uses in internal buffer for computing derivatives)
    
    Real evaluate( pimcConfigurations_t & configurations , std::array<int,2> timeSlices , int iChain );
    
    Real evaluate(pimcConfigurations_t & configurations); // evaluates the full action

    const auto & particleGroups () const {return _particleGroups;}

    private:
    geometryPBC_PIMC geo;
    int nChains;
    int nBeads;
    std::vector<particleGroup> _particleGroups ;
    Eigen::Tensor<Real, 3> distancesBuffer;
    Real tau;
    Real D;

};


template<class functor_t>
class potentialActionOneBody : public action
{
    public:

    potentialActionOneBody(Real tau_, functor_t V_ ,std::vector<particleGroup> groups_ ,geometryPBC_PIMC geo_): tau(tau_),V(V_),geo(geo_),_groups(groups_) {}

    Real evaluate(pimcConfigurations_t & configurations, std::array<int,2> timeRange, std::array<int,2>  particleRange  ) 
    {
        auto & data = configurations.dataTensor();

        auto sum = reduceOnPositions(V, data, timeRange, particleRange, configurations.getMask());

        return sum;
    }

    Real evaluate(pimcConfigurations_t & configurations, std::array<int,2> timeRange, int iChain  ) 
    {
        Real sum=0;
        for (const auto & group : _groups)
        {
            if ( (iChain >=group.iStart) and (iChain <= group.iStart) )
            {
                sum+=evaluate(configurations,timeRange,{iChain,iChain});
            }
        }
        return sum;
    }


    Real evaluate(pimcConfigurations_t & configurations ) 
    {
        const auto nChains = configurations.nChains();
        const auto nBeads = configurations.nBeads();
        Real sum=0;

        for (const auto & group : _groups)
        {
        
        sum += evaluate(configurations, {0,nBeads-1},{group.iStart,group.iEnd});
        }

        return sum;
    }

    private:
    functor_t V;
    Real tau;
    geometryPBC_PIMC geo;
    std::vector<particleGroup> _groups;

};






template<class functor_t>
class potentialActionTwoBody
{
    public:

    potentialActionTwoBody(Real tau_, int nChains_  , int nBeads_, functor_t V_, geometryPBC_PIMC geo_) : tau(tau_),V(V_),geo(geo_),nChains(nChains_),nBeads(nBeads_) ,
    bufferDistances(nChains_,getDimensions(),nBeads_)
     {}

     Real evaluate(const pimcConfigurations_t & configurations, std::array<int,2> timeRange, int iChain)
     {
         Real sum=0;

         for (const auto & [groupA, groupB] : _groups)
         {
             bool isInA = groupA.contains(iChain);
             bool isInB = groupB.contains(iChain);

            if ( isInA or isInB   )
            {
                auto & otherGroup = isInA ? groupB : groupA;

                sum+=evaluate(configurations,timeRange,iChain,{otherGroup.iStart, otherGroup.iEnd } );
            }


         }

         return sum;

     };

    Real evaluate(const pimcConfigurations_t & configurations, std::array<int,2> timeRange, int iChain, std::array<int,2>  particleRange  ) 
    {
        const auto & data = configurations.dataTensor();

        geo.updateEqualTimeDifferences(bufferDistances, data, timeRange , iChain, particleRange) ;

        auto sum=reduceOnDifferences(bufferDistances,timeRange,iChain,particleRange,configurations.getMask());
        


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


    Real reduceOnDifferences(const Eigen::Tensor<Real,3> & tn, std::array<int,2> timeRange, int j,  std::array<int,2> particleRange  ) const 
    {
         Real sum=0;

         auto iStartLeft = particleRange[0];
         auto iEndLeft = std::min(particleRange[1],j);

         auto iStartRight = iEndLeft + 1;
         auto iEndRight = particleRange[1];


        for(int t=timeRange[0];t<=timeRange[1] ; t++)
        {
            for (size_t i=iStartLeft;i<=iEndLeft;i++ )
                {
                    sum+=V( tn(i,0,t ) , tn(i,1,t) , tn(i,2,t)  ) ;
                }
        for (size_t i=iStartRight;i<=iEndRight;i++ )
                {
                    sum+=V( tn(i,0,t ) , tn(i,1,t) , tn(i,2,t)  ) ;
                }
        }

        return sum;

    }
    
    Real reduceOnDifferences(const Eigen::Tensor<Real,3> & tn, std::array<int,2> timeRange, int j,  std::array<int,2> particleRange  , const mask & mask) const 
    {
         Real sum=0;

         auto iStartLeft = particleRange[0];
         auto iEndLeft = std::min(particleRange[1],j);

         auto iStartRight = iEndLeft + 1;
         auto iEndRight = particleRange[1];


        for(int t=timeRange[0];t<=timeRange[1] ; t++)
        {
            for (size_t i=iStartLeft;i<=iEndLeft;i++ )
                {
                    sum+= ( mask(i,t) == 0 ? 0 : V( tn(i,0,t ) , tn(i,1,t) , tn(i,2,t)  ) ) ;
                }
        for (size_t i=iStartRight;i<=iEndRight;i++ )
                {
                    sum+=( mask(i,t)==0  ? 0 : V( tn(i,0,t ) , tn(i,1,t) , tn(i,2,t)  ) )  ;
                }
        }


        return sum;

    }
    



    functor_t V;
    Real tau;
    Eigen::Tensor<Real,3> bufferDistances;
    geometryPBC_PIMC geo;
    int nChains,nBeads;
    std::vector< std::array<particleGroup,2 >  > _groups;

};


class sumAction : public action
{
    public:
    using action::evaluate;



    sumAction( std::vector<action *> actions_) : _actions(actions_) {}

    virtual Real evaluate(pimcConfigurations_t & configurations, std::array<int,2> timeRange, int iChain)
    {
        Real sum=0;
        for ( auto S : _actions)
        {
            sum+=S->evaluate(configurations,timeRange , iChain);
        }
        return sum;
    }

    virtual Real evaluate( pimcConfigurations_t & pimcConfigurations)
    {
        Real sum=0;
        for(auto S : _actions)
        {
            sum+=S->evaluate(pimcConfigurations);
        }
        return sum;
    }

    auto operator[](int i) {return _actions[i];}




    auto  & getActions() {return _actions;}

    private:

    std::vector<action* > _actions;
};


class firstOrderAction : public sumAction
{
public:
    firstOrderAction( action * kineticAction, action * potentialAction) : sumAction::sumAction( {kineticAction,potentialAction}) 
    {

    }
    
    auto getKineticAction() {return (*this)[0];}
    auto getPotentialAction() {return (*this)[1];}

};


}

