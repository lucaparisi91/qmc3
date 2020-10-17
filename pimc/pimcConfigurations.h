
#ifndef PIMCCONFIGURATIONS
#define PIMCCONFIGURATIONS


#include <vector>
#include "../src/traits.h"
#include <unsupported/Eigen/CXX11/Tensor>

namespace pimc
{
    
class mask
{
public:
    mask( int nBeads, int nChains) : _mask(nChains,nBeads) {_mask.setConstant(1);}
    auto & operator()(int i,int j) {return _mask(i,j);}
    const auto &  operator()(int i,int j) const {return _mask(i,j);}

private:
    Eigen::Tensor<int, 2> _mask;
};


    struct particleGroup
    {
        bool contains(int iParticle) {return (iParticle>= iStart) and (iParticle<=iEnd);}
        int iStart;
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
        _data(nChains, dimensions,timeSlices), particleGroups(particleGroups_)   ,N(nChains),M(timeSlices),
        _mask(timeSlices,nChains)
        {
            
        };



        auto &  dataTensor() {return _data;}
         const auto &  dataTensor() const {return _data;}


        auto data() {return _data.data();}
        
        auto nChains() const {return N;} 
        
        auto nBeads() const {return M; } 

        const auto & getMask() const {return _mask;}

        void deleteBeads(   std::array<int,2> timeRange, int iChain );

        void createBeads(   std::array<int,2> timeRange, int iChain );


        private:

        std::vector<particleGroup> particleGroups;
        std::vector<worm> worms;
        configurationsStorage_t _data;
        int M;
        int N;
        mask _mask;

    };

    using configurations_t = pimcConfigurations; 

};

#endif