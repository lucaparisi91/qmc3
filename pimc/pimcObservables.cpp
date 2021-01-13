#include "pimcObservables.h"
namespace pimc
{

Real thermodynamicEnergyEstimator::operator()(configurations_t & confs, firstOrderAction & S)
{
    auto & geo = S.getGeometry();

    auto & kA = S.getKineticAction();
    auto & potA = S.getPotentialAction();
    
    auto sA=kA.evaluate(confs);
    auto sV=potA.evaluate(confs);

    auto beta = confs.nBeads() * kA.getTimeStep(); 
    sA/=beta*confs.nParticles();
    sV/=beta*confs.nParticles();

    return sV - sA +  getDimensions()/(2.*kA.getTimeStep());
}



Real virialEnergyEstimator::operator()(configurations_t & confs, firstOrderAction & S)
{
    auto & geo = S.getGeometry();
    auto & Spot = S.getPotentialAction();
    Real e=0;
    Real e2=0 , e3 = 0 , e4=0;

    int N=0;
    buffer.setConstant(0.);

    for ( const auto & group : confs.getGroups() )
    {
        N+=group.iEnd - group.iStart + 1;
        Spot.addGradient(confs,{0,confs.nBeads()-1},{group.iStart,group.iEnd},buffer);
    }
    // 
    const auto & data = confs.dataTensor(); 


    // compute rC
    rC.setConstant(0);
    for (const auto & group : confs.getGroups() )
    {
        

        for (int alpha=0; alpha< confs.nBeads();alpha++)
        {
            // sum on time slices greater than alpha, same chain
            for (int t=alpha;t<confs.nBeads();t++)
            {
                for (int i = group.iStart ; i<=group.iEnd ; i++)
                {
                    for (int d=0;d<getDimensions();d++)
                    {
                        rC(i,d,alpha)+=0.5*data(i,d,t)/ confs.nBeads();
                    }
                }
                
            }
            // sum on time slices greater than alpha, next chain
            for (int t=0;t<alpha;t++)
            {
                for (int i = group.iStart ; i<=group.iEnd ; i++)
                {
                    int ii=confs.getChain(i).next;

                    for (int d=0;d<getDimensions();d++)
                    {
                        rC(i,d,alpha)+=0.5*data(ii,d,t)/ confs.nBeads();
                    }
                }
            }

            // sum on time slices lesser than alpha, prev chain
            for (int t=alpha + 1;t< confs.nBeads();t++)
            {
                for (int i = group.iStart ; i<=group.iEnd ; i++)
                {
                    int ii=confs.getChain(i).prev;

                    for (int d=0;d<getDimensions();d++)
                    {
                        rC(i,d,alpha)+=0.5*data(ii,d,t)/ confs.nBeads();
                    }
                }
            }

            // sum on time slices lesser then alpha , same chain
            for (int t=0;t<=alpha;t++)
            {
                for (int i = group.iStart ; i<=group.iEnd ; i++)
                {
                    for (int d=0;d<getDimensions();d++)
                    {
                        rC(i,d,alpha)+=0.5*data(i,d,t)/ confs.nBeads();
                    }
                }
                
            }

           /*  std::cout << alpha<<" ";
                for (int i = group.iStart ; i<=group.iEnd ; i++)
                {
                    for (int d=0;d<getDimensions();d++)
                    {
                        std::cout << rC(i,d,alpha) << " ";
                    }
                }
            std::cout <<std::endl;
             */
            


        }

        // second term in the virial estimator
        // Does not give any contribution for a classical system


        for (int t=1;t<confs.nBeads();t++)
            {
                for (int i = group.iStart ; i<=group.iEnd ; i++)
                {
                    int inext = confs.getChain(i).next;

                    for(int d=0;d<getDimensions();d++)
                    {
                        e2+= 
                        /*
                        ( geo.difference( data(inext,d,t ) - data(i,d,t) ,d )  ) *
                        ( geo.difference( data(inext,d,t ) - data(i,d,t) ,d )  );
                        
                        geo.difference( -data(i,d,t) + data(inext,d,t ) ,d )*
                        geo.difference( -data(i,d,t+1) + data(i,d,t ) ,d );
                        */
                        (-data(i,d,t) + data(inext,d,t ) + data(i,d,confs.nBeads())  - data(inext,d,0)   )*
                        (  data(inext,d,t-1)   - data(inext,d,t)   );


                    }
                    
                    
                }
            }
        //std::cout << e2 << std::endl;

        {
            int t=0;
                for (int i = group.iStart ; i<=group.iEnd ; i++)
                {
                    int inext = confs.getChain(i).next;

                    for(int d=0;d<getDimensions();d++)
                    {
                        e2+= (-data(i,d,t) + data(i,d,confs.nBeads() )   )*( data(i,d,confs.nBeads()-1) - data(i,d,confs.nBeads()  ) );
                        
                    }
                    
                    
                }
            } 


        // third term in the virial estimator
        for (int t=0;t<confs.nBeads();t++)
            {
                for (int i = group.iStart ; i<=group.iEnd ; i++)
                {
                    for(int d=0;d<getDimensions();d++)
                    {
                        e3+=(data(i,d,t) - rC(i,d,t) )*buffer(i,d,t);
                    }
                }
             
            }
    }


    Real beta = S.getTimeStep() * confs.nBeads();
    e4= S.getPotentialAction().evaluate(confs);
    e4/=(beta *N);

    e3/=(2*N*beta);

    e2/=(2*N*beta*beta);
    
    // classical gas free contribution
    Real e1 = getDimensions()/(2*beta);
    return e1 + e2  + e3 + e4;

}



};