#ifndef TRAITS_H
#define TRAITS_H
#include <vector>
#include <unsupported/Eigen/CXX11/Tensor>


class geometry;

using real_t = double;
using geometry_t = geometry;
using state_t = Eigen::Tensor<real_t,2> ;
using states_t = std::vector<state_t> ;

using difference_t = state_t;
using differences_t = std::vector<difference_t>;

using distance_t = Eigen::Tensor<real_t,1> ;
using distances_t = std::vector<distance_t>;


#endif