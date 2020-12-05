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

    
    buffer.setConstant(0.);


    for ( const auto & group : confs.getGroups() )
    {
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
                for (int i = group.iStart ; i<group.iEnd ; i++)
                {
                    for (int d=0;d<getDimensions();d++)
                    {
                        rC(i,d,t)+=data(i,d,t);
                    }
                }
                
            }
            // sum on time slices greater than alpha, next chain
            for (int t=0;t<=alpha;t++)
            {
                for (int i = group.iStart ; i<=group.iEnd ; i++)
                {
                    int ii=confs.getChain(i).next;

                    for (int d=0;d<getDimensions();d++)
                    {
                        rC(i,d,t)+=data(ii,d,t);
                    }
                }
            }

            // sum on time slices lesset than alpha, prev chain
            for (int t=confs.nBeads() -1 - alpha;t< confs.nBeads();t++)
            {
                for (int i = group.iStart ; i<=group.iEnd ; i++)
                {
                    int ii=confs.getChain(i).prev;

                    for (int d=0;d<getDimensions();d++)
                    {
                        rC(i,d,t)+=data(ii,d,t);
                    }
                }
            }

            // sum on time slices lesser then alpha , same chain
            for (int t=0;t<=alpha;t++)
            {
                for (int i = group.iStart ; i<group.iEnd ; i++)
                {
                    for (int d=0;d<getDimensions();d++)
                    {
                        rC(i,d,t)+=data(i,d,t);
                    }
                }
                
            }



        }

    }

    // sum on neighbouring particle differences




}



};