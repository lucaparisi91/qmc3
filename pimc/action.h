#ifndef ACTION_H
#define ACTION_H

#include "unsupported/Eigen/CXX11/Tensor"
#include "geometryPMC.h"
#include "pimcConfigurations.h"
#include "tools.h"
#include "toolsPimc.h"
#include "qmcExceptions.h"


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
                    sum+= ( mask(t,i)==0 ? 0 :

                    #if DIMENSIONS == 3
                     V( tn( t,0, i  ) , tn(t,1,i) , tn(t,2,i)  ) 
                    #endif

                    #if DIMENSIONS == 1
                     V( tn( t,0, i  )   ) 
                    #endif

                    #if DIMENSIONS == 2
                     V( tn( t,0, i  ) , tn(t,1,i)   ) 
                    #endif

                     );
                }
        return sum;
    };



/*
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
*/

    template<class functor_t>
    Real reduceOnPositions(const functor_t & V,const Eigen::Tensor<Real,3> & tn, std::array<int ,2 > timeRange, std::array<int, 2>  particleRange, const mask & mask)
    {
        Real sum=0;

         
            for(int t=timeRange[0]+1;t<=timeRange[1]-1 ; t++ )
            {
              for (size_t i=particleRange[0];i<=particleRange[1];i++ )
                {
                    sum+= 
                    #if DIMENSIONS == 3
                     V( tn( i,0, t  ) , tn(i,1,t) , tn(i,2,t) )
                     #endif
                     #if DIMENSIONS == 1
                     V( tn( i,0, t  )) 
                     #endif
                     #if DIMENSIONS == 2
                     V( tn( i,0, t  ),tn( i,1, t  )) 
                     #endif

                       ;
                }
            }


            if (timeRange[1] > timeRange[0] )
            {
                int t = timeRange[0];
                for (size_t i=particleRange[0];i<=particleRange[1];i++ )
                {
                    sum+= 
                    #if DIMENSIONS == 3
                     0.5*V( tn( i,0, t  ) , tn(i,1,t) , tn(i,2,t) )
                     #endif
                     #if DIMENSIONS == 1
                     0.5*V( tn( i,0, t  )) 
                     #endif
                     #if DIMENSIONS == 2
                     0.5*V( tn( i,0, t  ),tn( i,1, t  )) 
                     #endif

                      
                     ;
                }

            }



            if (timeRange[1] >= timeRange[0] )
            {
            

                int t = timeRange[1];
                for (size_t i=particleRange[0];i<=particleRange[1];i++ )
                {
                    sum+=  
                    #if DIMENSIONS == 3
                     0.5*V( tn( i,0, t  ) , tn(i,1,t) , tn(i,2,t) )
                     #endif
                     #if DIMENSIONS == 1
                     0.5*V( tn( i,0, t  )) 
                     #endif
                     #if DIMENSIONS == 2
                     0.5*V( tn( i,0, t  ),tn( i,1, t  )) 
                     #endif

                      
                     ;
                }

                t = timeRange[1]+1;
                for (size_t i=particleRange[0];i<=particleRange[1];i++ )
                {
                    sum+= 
                    #if DIMENSIONS == 3
                     0.5*V( tn( i,0, t  ) , tn(i,1,t) , tn(i,2,t) )
                     #endif
                     #if DIMENSIONS == 1
                     0.5*V( tn( i,0, t  )) 
                     #endif
                     #if DIMENSIONS == 2
                     0.5*V( tn( i,0, t  ),tn( i,1, t  )) 
                     #endif

                      
                     ;
                }
            }               
                   


        return sum;
    };

    
using pimcConfigurations_t = pimcConfigurations;

class action
{
    public:


    action(Real timeStep, const geometryPBC_PIMC & geo_) : _geo(geo_), _timeStep(timeStep) {}


    virtual Real evaluate(configurations_t & configurations, std::array<int,2> timeRange, int iParticle)=0;

    virtual Real evaluate(pimcConfigurations_t & configurations, std::array<int,2> timeRange, int iChain1 , int iChain2  )=0; 

    virtual Real evaluate( configurations_t & pimcConfigurations)=0; // evaluates the whole action


    virtual void addGradient(const configurations_t & pimcConfigurations,const std::array<int,2> & timeRange,const  std::array<int,2> & particleRange, const Eigen::Tensor<Real,3> & gradientBuffer){throw missingImplementation("Gradient not implemented for this action.");}





    const auto & getGeometry() const {return _geo ;}
    auto & getGeometry() {return _geo ;}

    auto getTimeStep() const {return _timeStep;}

    private:
    geometryPBC_PIMC _geo;
    Real _timeStep;

};

class kineticAction : public action
{
    public:
    using action::evaluate;

    kineticAction(Real tau_, int nChains_ , int nBeads_, const geometryPBC_PIMC & geo_
    ) : nChains(nChains_),nBeads(nBeads_),distancesBuffer(nBeads_, getDimensions( ) , nChains_),
    D(0.5)  , action::action(tau_,geo_) {}

    virtual Real evaluate( pimcConfigurations_t & configurations , std::array<int,2> timeSlices , std::array<int,2> chainRange ); // evaluates the kinetic action to recalculate from beads(iParticle, timeSliceStart: timeSliceEnd, uses in internal buffer for computing derivatives)
    
    virtual Real evaluate( pimcConfigurations_t & configurations , std::array<int,2> timeSlices , int iChain );

    virtual Real evaluate( pimcConfigurations_t & configurations , std::array<int,2> timeSlices , int iChain1 , int iChain2 );
    
    Real evaluate(pimcConfigurations_t & configurations); // evaluates the full action



    


    private:
    int nChains;
    int nBeads;

    Eigen::Tensor<Real, 3> distancesBuffer;
    Real D;

};

template<class V_t,class gradX_t >
class potentialFunctor{
public:
    potentialFunctor(V_t V_,gradX_t gradX_) : V(V_),_gradX(gradX_) {}

    Real operator()(Real x) const {return V(x);}
    Real gradX(Real x) const  {return _gradX(x);}

private:
    V_t V;
    gradX_t _gradX;
};

template<class V_t,class gradX_t >
auto makePotentialFunctor(V_t V_,gradX_t X_)
{
    return potentialFunctor<V_t,gradX_t>(V_,X_);    
}

template<class functor_t>
class potentialActionOneBody : public action
{
    public:

    potentialActionOneBody(Real tau_, functor_t V_ ,geometryPBC_PIMC geo_): V(V_),action::action(tau_,geo_) {}

    Real evaluate(pimcConfigurations_t & configurations, std::array<int,2> timeRange, std::array<int,2>  particleRange  ) 
    {
        auto & data = configurations.dataTensor();

        auto sum = reduceOnPositions(V, data, timeRange, particleRange, configurations.getMask());

        return getTimeStep()*sum;

    }

    Real evaluate(pimcConfigurations_t & configurations, std::array<int,2> timeRange, int iChain  ) 
    {
        Real sum=0;
        const auto & groups = configurations.getGroups();
        for (const auto & group : groups)
        {
            if ( (iChain >=group.iStart) and (iChain <= group.iEnd) )
            {
                sum+=evaluate(configurations,timeRange,{iChain,iChain});
            }
        }
        return sum;
    }

    Real evaluate(pimcConfigurations_t & configurations, std::array<int,2> timeRange, int iChain1 , int iChain2  ) 
    {
        return evaluate(configurations,timeRange,iChain1) + evaluate(configurations,timeRange,iChain2); 

    }

    Real evaluate(pimcConfigurations_t & configurations ) 
    {
        const auto nChains = configurations.nChains();
        const auto nBeads = configurations.nBeads();
        Real sum=0;
        const auto & groups = configurations.getGroups();
        for (const auto & group : groups)
        {
        
        sum += evaluate(configurations, {0,nBeads-1},{group.iStart,group.iEnd});
        }

        return sum;
    }


     virtual void addGradient(const configurations_t & pimcConfigurations,const std::array<int,2> & timeRange, const std::array<int,2> & particleRange,  Eigen::Tensor<Real,3> & gradientBuffer){

         const auto & data = pimcConfigurations.dataTensor();

         for (int t=timeRange[0];t<=timeRange[1];t++)
            for (int i=particleRange[0] ; i<=particleRange[1];i++ )
         {
             #if DIMENSIONS == 1
            gradientBuffer(i,0,t)+=V.gradX( data(i,0,t)  );
             #endif

             #if DIMENSIONS == 2
            gradientBuffer(i,0,t)+=gradV.X( data(i,0,t) ,  data(i,1,t)  );
            gradientBuffer(i,1,t)+=gradV.Y( data(i,0,t) ,  data(i,1,t)  );
            
             #endif


            #if DIMENSIONS == 3
            gradientBuffer(i,0,t)+=gradV.X( data(i,0,t) ,  data(i,1,t) , data(i,2,t) );
            gradientBuffer(i,1,t)+=gradV.Y( data(i,0,t) ,  data(i,1,t) , data(i,2,t) );
            gradientBuffer(i,2,t)+=gradV.Z( data(i,0,t) ,  data(i,1,t) , data(i,2,t) );            
             #endif

         }
     }



    private:
    functor_t V;
};


template<class functor_t>
class potentialActionTwoBody
{
    public:

    potentialActionTwoBody(Real tau_, int nChains_  , int nBeads_, functor_t V_, geometryPBC_PIMC geo_) : tau(tau_),V(V_),geo(geo_),nChains(nChains_),nBeads(nBeads_) ,
    bufferDistances(nChains_,getDimensions(),nBeads_) 
     {}

     Real evaluate( pimcConfigurations_t & configurations, std::array<int,2> timeRange, int iChain)
     {
         Real sum=0;

         auto & [groupA, groupB] = _groups[0];

         
        bool isInA = groupA.contains(iChain);
        bool isInB = groupB.contains(iChain);

        if ( isInA or isInB   )
        {
            auto & otherGroup = isInA ? groupB : groupA;

            sum+=evaluate(configurations,timeRange,iChain,{otherGroup.iStart, otherGroup.iEnd } );
        }         

         return sum;

     };


     Real evaluate(pimcConfigurations_t & configurations, std::array<int,2> timeRange, int iChain1, int iChain2)
     {
        auto & [groupA, groupB] = _groups[0];

         // evaluates were all particles are
        bool is1inA = groupA.contains(iChain1);
        bool is1inB = groupB.contains(iChain1);
        bool is2inA = groupA.contains(iChain2);
        bool is2inB = groupB.contains(iChain2);

        bool is1inAorB =  is1inA or is1inB;
        bool is2inAorB = is2inA or is2inB;
        

         Real sum=0;

        // both particles are contained in setA or setB
         if ( is1inAorB and is2inAorB   )
         {
            const auto & groupLeft =  is1inA ? groupB : groupA;
            const auto & groupRight = is2inA ? groupB : groupA;

            return evaluate(configurations,timeRange,iChain1,iChain2, {groupLeft.iStart,groupLeft.iEnd},{groupRight.iStart,groupRight.iEnd}     );
        }
        else
        {
            auto iChain = is1inAorB ? iChain1 : iChain2; 
            return evaluate(configurations,timeRange,iChain);
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

    Real evaluate(const pimcConfigurations_t & configurations, std::array<int,2> timeRange, int iChain1, int iChain2, std::array<int,2>  particleRange1, std::array<int,2> particleRange2  ) 
    {
        auto sum1=evaluate(configurations,timeRange,iChain1,particleRange1);
        // sum distances changed with particle j . Avoid multiple-counting of particle (iChain1,iChain2)

        std::array<int,2> firstInterval = {particleRange2[0],std::min(particleRange2[1],iChain1)};

        std::array<int,2> secondInterval = {firstInterval[1]+1,particleRange2[1]};

        auto sum2Left=evaluate(configurations,timeRange,iChain2,firstInterval);
        auto sum2Right=evaluate(configurations,timeRange,iChain2,secondInterval);


        return sum1 + sum2Left + sum2Right;
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
                    sum+= ( mask(t,i) == 0 ? 0 : V( tn(i,0,t ) , tn(i,1,t) , tn(i,2,t)  ) ) ;
                }
        for (size_t i=iStartRight;i<=iEndRight;i++ )
                {
                    sum+=( mask(t,i)==0  ? 0 : V( tn(i,0,t ) , tn(i,1,t) , tn(i,2,t)  ) )  ;
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


    sumAction(  std::vector<action *> actions_) : _actions(actions_),action::action( actions_[0]->getTimeStep() ,    actions_[0]->getGeometry()  ) {}

    virtual Real evaluate(pimcConfigurations_t & configurations, std::array<int,2> timeRange, int iChain)
    {
        Real sum=0;
        for ( auto S : _actions)
        {
            sum+=S->evaluate(configurations,timeRange , iChain);
        }
        return sum;
    }

    virtual Real evaluate(pimcConfigurations_t & configurations, std::array<int,2> timeRange, int iChain1 , int iChain2)
    {
        Real sum=0;
        for ( auto S : _actions)
        {
            sum+=S->evaluate(configurations,timeRange , iChain1, iChain2);
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

    auto & operator[](int i) {return *(_actions[i]);}

    auto  & getActions() {return _actions;}

    private:

    std::vector<action* > _actions;
};

class firstOrderAction : public sumAction
{
public:
    firstOrderAction(  action * kineticAction, action * potentialAction) : sumAction::sumAction( {kineticAction,potentialAction}) 
    {

    }

    auto & getKineticAction() {return (*this)[0];}
    auto & getPotentialAction() {return (*this)[1];}

};


}


#endif