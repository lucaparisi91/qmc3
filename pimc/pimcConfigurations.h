
#ifndef PIMCCONFIGURATIONS
#define PIMCCONFIGURATIONS


#include <vector>
#include "../src/traits.h"
#include <unsupported/Eigen/CXX11/Tensor>

namespace pimc
{
    
    struct particleGroup
    {
        int istart;
        int iEnd;
        Real mass;
    };

    struct worm
    {
        int iChain;
        int iTail;
        int iHead;
    };

    class pimcConfigurations
    {
        public:
        using configurationsStorage_t =  Eigen::Tensor<Real, 3> ;

        pimcConfigurations(size_t timeSlices, size_t nChains, int dimensions, std::vector<particleGroup>  particleGroups_) : 
        _data(nChains, dimensions,timeSlices), particleGroups(particleGroups_)   ,N(nChains),M(timeSlices)
        {
            
        };

        auto &  dataTensor() {return _data;}
         const auto &  dataTensor() const {return _data;}


        auto data() {return _data.data();}
        
        auto nChains() const {return N;} 
        
        auto nBeads() const {return M; } 

        private:

        std::vector<particleGroup> particleGroups;
        std::vector<worm> worms;
        configurationsStorage_t _data;
        int M;
        int N;


    };


};

#endif