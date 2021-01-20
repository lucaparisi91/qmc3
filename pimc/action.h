#ifndef ACTION_H
#define ACTION_H

#include "unsupported/Eigen/CXX11/Tensor"
#include "geometryPMC.h"
#include "pimcConfigurations.h"
#include "tools.h"
#include "toolsPimc.h"
#include "qmcExceptions.h"
#include <memory>


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
                    sum+=

                    #if DIMENSIONS == 3
                     V( tn( t,0, i  ) , tn(t,1,i) , tn(t,2,i)  ) 
                    #endif

                    #if DIMENSIONS == 1
                     V( tn( t,0, i  )   ) 
                    #endif

                    #if DIMENSIONS == 2
                     V( tn( t,0, i  ) , tn(t,1,i)   ) 
                    #endif

                     ;
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

        for(int t=timeRange[0]+1;t<=timeRange[1] ; t++ )
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



            
            {

                int t = timeRange[1]+1;
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

    action(){}

    action(Real timeStep, const geometryPBC_PIMC & geo_) : _geo(geo_), _timeStep(timeStep) {}


    virtual Real evaluate(configurations_t & configurations, std::array<int,2> timeRange, int iParticle)=0;

    virtual Real evaluate( configurations_t & pimcConfigurations)=0; // evaluates the whole action


    virtual void addGradient(const configurations_t & pimcConfigurations,const std::array<int,2> & timeRange,const  std::array<int,2> & particleRange,  Eigen::Tensor<Real,3> & gradientBuffer){throw missingImplementation("Gradient not implemented for this action.");}



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


#if DIMENSIONS == 1
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

#endif

#if DIMENSIONS == 3
template<class V_t,class gradX_t , class gradY_t,class gradZ_t>
class potentialFunctor{
public:
    
    potentialFunctor(V_t V_,gradX_t gradX_,gradY_t gradY_,gradZ_t gradZ_) : V(V_),_gradX(gradX_),_gradY(gradY_),_gradZ(gradZ_) {}

    Real operator()(Real x,Real y , Real z) const { return V(x,y,z); }

    
    Real gradX(Real x,Real y,Real z) const  {return _gradX(x,y,z);}
    Real gradY(Real x,Real y,Real z) const  {return _gradY(x,y,z);}
    Real gradZ(Real x,Real y,Real z) const  {return _gradZ(x,y,z);}



private:
    V_t V;
    gradX_t _gradX;
    gradY_t _gradY;
    gradZ_t _gradZ;
};


template<class V_t,class gradX_t ,class gradY_t , class gradZ_t>
auto makePotentialFunctor(V_t V_,gradX_t X_,gradY_t Y_, gradZ_t Z_)
{
    return potentialFunctor<V_t,gradX_t,gradY_t,gradZ_t>(V_,X_,Y_,Z_);
}

#endif



template<class functor_t>
class potentialActionOneBody : public action
{
    public:
    potentialActionOneBody(Real tau, geometryPBC_PIMC geo_, const json_t & j):
    potentialActionOneBody::potentialActionOneBody(tau,functor_t(j["potential"]) , geo_   ) {} 


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

    virtual void addGradient(const configurations_t & pimcConfigurations,const std::array<int,2> & timeRange,const  std::array<int,2> & particleRange,  Eigen::Tensor<Real,3> & gradientBuffer)
   {

         const auto & data = pimcConfigurations.dataTensor();

         for (int t=timeRange[0];t<=timeRange[1];t++)
            for (int i=particleRange[0] ; i<=particleRange[1];i++ )
         {
             #if DIMENSIONS == 1
            gradientBuffer(i,0,t)+=V.gradX( data(i,0,t)  )*getTimeStep();
             #endif

             #if DIMENSIONS == 2
            gradientBuffer(i,0,t)+=gradV.X( data(i,0,t) ,  data(i,1,t)  )*getTimeStep();
            gradientBuffer(i,1,t)+=gradV.Y( data(i,0,t) ,  data(i,1,t)  )*getTimeStep();
            
             #endif


            #if DIMENSIONS == 3
            gradientBuffer(i,0,t)+=V.gradX( data(i,0,t) ,  data(i,1,t) , data(i,2,t) )*getTimeStep();
            gradientBuffer(i,1,t)+=V.gradY( data(i,0,t) ,  data(i,1,t) , data(i,2,t) )*getTimeStep();
            gradientBuffer(i,2,t)+=V.gradZ( data(i,0,t) ,  data(i,1,t) , data(i,2,t) )*getTimeStep();            
             #endif
         }

       /*  // t[0] only counts half
        {
            int t=timeRange[0];
            for (int i=particleRange[0] ; i<=particleRange[1];i++ )
            {
                #if DIMENSIONS == 1
                gradientBuffer(i,0,t)+=V.gradX( data(i,0,t)  )*getTimeStep();
                #endif

                #if DIMENSIONS == 2
                gradientBuffer(i,0,t)+=0.5*gradV.X( data(i,0,t) ,  data(i,1,t)  )*getTimeStep();
                gradientBuffer(i,1,t)+=0.5*gradV.Y( data(i,0,t) ,  data(i,1,t)  )*getTimeStep();                
                #endif

                #if DIMENSIONS == 3
                gradientBuffer(i,0,t)+=0.5*V.gradX( data(i,0,t) ,  data(i,1,t) , data(i,2,t) )*getTimeStep();
                gradientBuffer(i,1,t)+=0.5*V.gradY( data(i,0,t) ,  data(i,1,t) , data(i,2,t) )*getTimeStep();
                gradientBuffer(i,2,t)+=0.5*V.gradZ( data(i,0,t) ,  data(i,1,t) , data(i,2,t) )*getTimeStep();            
                #endif
            }
        }
 */
      /*   // t[1] + 1 counts half
        {
            int t=timeRange[1] + 1;

            for (int i=particleRange[0] ; i<=particleRange[1];i++ )
            {
                #if DIMENSIONS == 1
                gradientBuffer(i,0,t)+=0.5*V.gradX( data(i,0,t)  )*getTimeStep();
                #endif

                #if DIMENSIONS == 2
                gradientBuffer(i,0,t)+=0.5*gradV.X( data(i,0,t) ,  data(i,1,t)  )*getTimeStep();
                gradientBuffer(i,1,t)+=0.5*gradV.Y( data(i,0,t) ,  data(i,1,t)  )*getTimeStep();                
                #endif

                #if DIMENSIONS == 3
                gradientBuffer(i,0,t)+=0.5*V.gradX( data(i,0,t) ,  data(i,1,t) , data(i,2,t) )*getTimeStep();
                gradientBuffer(i,1,t)+=0.5*V.gradY( data(i,0,t) ,  data(i,1,t) , data(i,2,t) )*getTimeStep();
                gradientBuffer(i,2,t)+=0.5*V.gradZ( data(i,0,t) ,  data(i,1,t) , data(i,2,t) )*getTimeStep();            
                #endif
            }
        } */
     }

    private:
    functor_t V;
};

template<class functor_t>
class potentialActionTwoBody : public action
{
    public:

    potentialActionTwoBody(Real tau_, int nChains_  , int nBeads_, geometryPBC_PIMC geo_, const json_t & j) :
    potentialActionTwoBody(tau_,nChains_,nBeads_,functor_t(j["potential"]),geo_,j["groupA"].get<int>() , j["groupB"].get<int>()  ) {}


    potentialActionTwoBody(Real tau_, int nChains_  , int nBeads_, functor_t V_, geometryPBC_PIMC geo_, int  iParticleGroupA_, int iParticleGroupB_) : V(V_),nChains(nChains_),nBeads(nBeads_) ,
    bufferDistances(nChains_,getDimensions(),nBeads_+1),
    iParticleGroupA(iParticleGroupA_),iParticleGroupB(iParticleGroupB_),
    action::action(tau_,geo_)
     {}

    virtual Real evaluate(configurations_t & configurations, std::array<int,2> timeRange, int iChain)
     {
         Real sum=0;


        const auto & groupA = configurations.getGroups()[iParticleGroupA];
        const auto & groupB = configurations.getGroups()[iParticleGroupB];

         
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
        throw missingImplementation("Two particles updates are not supported");
     };

    Real evaluate(const pimcConfigurations_t & configurations, std::array<int,2> timeRange, int iChain, std::array<int,2>  particleRange  ) 
    {
        const auto & geo = getGeometry();
        Real sum=0;

        const auto & data = configurations.dataTensor();

        for(int t=timeRange[0]+1;t<=timeRange[1];t++)
        {
            for(int j=particleRange[0];j<=particleRange[1];j++)
            {
                for(int d=0;d<getDimensions();d++)
                {
                    bufferDistances(j,d,t)=geo.difference( 
                        data(iChain,d,t) - data(j,d,t) ,d
                    );
                }

                #if DIMENSIONS == 3
                sum+=V(bufferDistances(j,0,t) , bufferDistances(j,1,t) , bufferDistances(j,2   ,t)  );
                #endif 

                 #if DIMENSIONS == 2
                sum+=V(bufferDistances(j,0,t) , bufferDistances(j,1,t)   );
                #endif 

                 #if DIMENSIONS == 1
                sum+=V(bufferDistances(j,0,t)  );
                #endif 

            }
        }

        {
            int t=timeRange[0];
            for(int j=particleRange[0];j<=particleRange[1];j++)
            {
                for(int d=0;d<getDimensions();d++)
                {
                    bufferDistances(j,d,t)=geo.difference( 
                        data(iChain,d,t) - data(j,d,t) ,d
                    );
                }

                #if DIMENSIONS == 3
                sum+=0.5*V(bufferDistances(j,0,t) , bufferDistances(j,1,t) , bufferDistances(j,2   ,t)  );
                #endif 

                 #if DIMENSIONS == 2
                sum+=0.5*V(bufferDistances(j,0,t) , bufferDistances(j,1,t)   );
                #endif 

                 #if DIMENSIONS == 1
                sum+=0.5*V(bufferDistances(j,0,t)  );
                #endif 

            }
        }


        {
            int t=timeRange[1]+1;
            for(int j=particleRange[0];j<=particleRange[1];j++)
            {
                for(int d=0;d<getDimensions();d++)
                {
                    bufferDistances(j,d,t)=geo.difference( 
                        data(iChain,d,t) - data(j,d,t) ,d
                    );
                }

                #if DIMENSIONS == 3
                sum+=0.5*V(bufferDistances(j,0,t) , bufferDistances(j,1,t) , bufferDistances(j,2   ,t)  );
                #endif 

                 #if DIMENSIONS == 2
                sum+=0.5*V(bufferDistances(j,0,t) , bufferDistances(j,1,t)   );
                #endif 

                 #if DIMENSIONS == 1
                sum+=0.5*V(bufferDistances(j,0,t)  );
                #endif 

            }
        }






        return sum*getTimeStep();
    }

    Real evaluate(pimcConfigurations_t & configurations)
    {
        if (iParticleGroupA == iParticleGroupB)
        {
            return evaluateOnSameGroup(configurations);
        }
        else
        {
            return  evaluateOnDifferentGroup(configurations);
        }

    };


    virtual void addGradient(const configurations_t & pimcConfigurations,const std::array<int,2> & timeRange,const  std::array<int,2> & particleRange,  Eigen::Tensor<Real,3> & gradientBuffer)
    {
        const auto & geo = getGeometry();

        // should only used in the diagonal sector with time periodic boundary conditions

         const auto & data = pimcConfigurations.dataTensor();
         for (int t=timeRange[0];t<=timeRange[1];t++)
            for (int i=particleRange[0] ; i<=particleRange[1];i++ )
                for (int j=particleRange[0] ; j<i;j++ )
                {

                    for(int d=0;d<getDimensions();d++)
                    {
                        bufferDistances(j,d,t)=geo.difference( 
                            data(i,d,t) - data(j,d,t) ,d
                        );
                    }

                    #if DIMENSIONS == 1
                    Real tmp=V.gradX( bufferDistances(j,0,t)  )*getTimeStep();

                    gradientBuffer(i,0,t)+=tmp;
                    gradientBuffer(j,0,t)-=tmp;
                    #endif

                    #if DIMENSIONS == 3
                    Real tmp0=V.gradX( bufferDistances(j,0,t), bufferDistances(j,1,t), bufferDistances(j,2,t)  )*getTimeStep();
                    Real tmp1=V.gradY( bufferDistances(j,0,t), bufferDistances(j,1,t), bufferDistances(j,2,t)  )*getTimeStep();
                    Real tmp2=V.gradZ( bufferDistances(j,0,t), bufferDistances(j,1,t), bufferDistances(j,2,t)  )*getTimeStep();


                    gradientBuffer(i,0,t)+=tmp0;
                    gradientBuffer(j,0,t)-=tmp0;
                    gradientBuffer(i,1,t)+=tmp1;
                    gradientBuffer(j,1,t)-=tmp1;
                    gradientBuffer(i,2,t)+=tmp2;
                    gradientBuffer(j,2,t)-=tmp2;

                    #endif
                }
    


  

   

    }


    private:


     Real evaluateOnDifferentGroup(const pimcConfigurations_t & configurations)
     {
         throw missingImplementation("evaluateOnDifferentGroup");

     }

    Real evaluateOnSameGroup(const pimcConfigurations_t & configurations)
    {
        Real sum=0;

        const auto & geo = getGeometry();

        const auto & data = configurations.dataTensor();


        const auto & group = configurations.getGroups()[iParticleGroupA];

            for(int t=1;t<configurations.nBeads();t++)
            {

            for(int i=group.iStart; i<=group.iEnd;i++)
            {
                for(int j=group.iStart;j< i ;j++)
                {
                    for(int d=0;d<getDimensions();d++)
                    {
                        bufferDistances(j,d,t)=geo.difference( 
                            data(i,d,t) - data(j,d,t) ,d
                        );
                    }

                    #if DIMENSIONS == 3
                    sum+=V(bufferDistances(j,0,t) , bufferDistances(j,1,t) , bufferDistances(j,2   ,t)  );
                    #endif 

                    #if DIMENSIONS == 2
                    sum+=V(bufferDistances(j,0,t) , bufferDistances(j,1,t)   );
                    #endif 

                    #if DIMENSIONS == 1
                    sum+=V(bufferDistances(j,0,t)  );
                    #endif 

                }

            }
         }


            {
                int t=0;
            for(int i=group.iStart; i<=group.iEnd;i++)
            {
                for(int j=group.iStart;j< i ;j++)
                {
                    for(int d=0;d<getDimensions();d++)
                    {
                        bufferDistances(j,d,t)=geo.difference( 
                            data(i,d,t) - data(j,d,t) ,d
                        );
                    }

                    #if DIMENSIONS == 3
                    sum+=0.5*V(bufferDistances(j,0,t) , bufferDistances(j,1,t) , bufferDistances(j,2   ,t)  );
                    #endif 

                    #if DIMENSIONS == 2
                    sum+=V(bufferDistances(j,0,t) , bufferDistances(j,1,t)   );
                    #endif 

                    #if DIMENSIONS == 1
                    sum+=0.5*V(bufferDistances(j,0,t)  );
                    #endif 

                }

            }
         }

        {
            int t=configurations.nBeads();
            for(int i=group.iStart; i<=group.iEnd;i++)
            {
                for(int j=group.iStart;j< i ;j++)
                {
                    for(int d=0;d<getDimensions();d++)
                    {
                        bufferDistances(j,d,t)=geo.difference( 
                            data(i,d,t) - data(j,d,t) ,d
                        );
                    }

                    #if DIMENSIONS == 3
                    sum+=0.5*V(bufferDistances(j,0,t) , bufferDistances(j,1,t) , bufferDistances(j,2   ,t)  );
                    #endif 

                    #if DIMENSIONS == 2
                    sum+=V(bufferDistances(j,0,t) , bufferDistances(j,1,t)   );
                    #endif 

                    #if DIMENSIONS == 1
                    sum+=0.5*V(bufferDistances(j,0,t)  );
                    #endif 

                }

            }
        } 



        return sum*getTimeStep();
    };


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
                    sum+= V( tn(i,0,t ) , tn(i,1,t) , tn(i,2,t)   ) ;
                }
        for (size_t i=iStartRight;i<=iEndRight;i++ )
                {
                    sum+= V( tn(i,0,t ) , tn(i,1,t) , tn(i,2,t)   )  ;
                }
        }


        return sum;

    }


    

    functor_t V;
   
    Eigen::Tensor<Real,3> bufferDistances;
    int nChains,nBeads;


    int iParticleGroupA;
    int iParticleGroupB;

};


class sumAction : public action
{
    public:
    using action::evaluate;

    sumAction(){}

    sumAction(  std::vector<std::shared_ptr<action> > actions_) : _actions(actions_),action::action( actions_[0]->getTimeStep() ,    actions_[0]->getGeometry()  ) {}

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

    virtual void addGradient(const configurations_t & pimcConfigurations,const std::array<int,2> & timeRange,const  std::array<int,2> & particleRange,  Eigen::Tensor<Real,3> & gradientBuffer)
    {
        for(auto S : _actions)
        {
            S->addGradient(pimcConfigurations,timeRange,particleRange,gradientBuffer);

        }
        

    }

    auto & operator[](int i) {return *(_actions[i]);}

    auto  & getActions() {return _actions;}

    private:
    
    std::vector<std::shared_ptr<action> > _actions;
};

class firstOrderAction : public sumAction
{
public:
    firstOrderAction() {}

    
    firstOrderAction(  std::shared_ptr<action> sKinetic, std::shared_ptr<action> sPot) : sumAction::sumAction( {sKinetic,sPot}) 
    {

    }

    auto & getKineticAction() {return (*this)[0];}
    auto & getPotentialAction() {return (*this)[1];}

};




/* Action constructor */



template<class T>
std::shared_ptr<pimc::action> __createAction(Real timeStep,geometryPBC_PIMC & geo,int nChains, int nBeads, const json_t & j) {
    return std::make_shared<T >(timeStep, geo,j);
    }

template<class T>
std::shared_ptr<pimc::action> __createActionWithBuffers(Real timeStep,geometryPBC_PIMC & geo, int nChains , int nBeads , const json_t & j) {
    return std::make_shared<T >(timeStep, nChains, nBeads, geo,j);
    }




class actionConstructor
{
    public:
    using geometry_t=geometryPBC_PIMC;


    actionConstructor(geometry_t geo, Real timeStep, int nChains_, int nBeads_) : _geo(geo),_timeStep(timeStep) , 
    _nChains(nChains_), _nBeads(nBeads_)
    {}


    auto createActionKey(const json_t & j)
    {
        return j["kind"].get<std::string>() + "_" + j["potential"]["kind"].get<std::string>();


    }

    std::shared_ptr<pimc::action> createAction(const json_t & j)
    {
        // creates a copy of the json input and supplement with additional essential information
     

        std::string key=createActionKey(j) ;


        creatorMap_t::const_iterator i;
        i=creatorMap.find(key);
        if (i!= creatorMap.end())
        {
	    return (i->second)( _timeStep,_geo,_nChains,_nBeads,j  );
        }
         else
        {
	    throw factoryIdNotRecorded(key);
        }
    }


    std::vector<std::shared_ptr<pimc::action> > createActions(const json_t & j)
    {
        std::vector< std::shared_ptr<pimc::action> > actions;
        for (const auto & jAction : j)
        {
            actions.push_back(createAction(jAction));
        }

        return actions;

    }





    template<class V_t>
    void registerPotential()
    {
        
        std::string key= "oneBody_" + V_t::name() ;
        creatorMap[key]= &__createAction<potentialActionOneBody<V_t> > ;

        key= "twoBody_" + V_t::name() ;
        creatorMap[key]= &__createActionWithBuffers<potentialActionTwoBody<V_t> > ;



    }

    

    private:

    geometry_t _geo;
    Real _timeStep;

    int _nChains;
    int _nBeads;

    

    typedef std::shared_ptr<pimc::action> (*actionCreatorFunc) (Real timeStep,geometryPBC_PIMC & geo, int nChains,int nBeads,const json_t & j);
    using creatorMap_t = std::map<std::string,actionCreatorFunc>;
    creatorMap_t creatorMap;
};










}






#endif