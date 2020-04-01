#ifndef TABLE_DISTANCES_H
#define TABLE_DISTANCES_H

#include "traits.h"
#include <unsupported/Eigen/CXX11/Tensor>
#include <map>
#include <unordered_map>

class tableDistances
{
public:
	using diff_t = Eigen::Tensor<real_t,2>;
	using distances_t = Eigen::Tensor<real_t ,1>;
	using state_t = Eigen::Tensor<real_t,2>;
	using states_t = std::vector<state_t>;
	tableDistances(geometry_t & geo_) : geo(&geo_){}

	void update(const states_t & states); // update all registered distances

	void add(int setA );  
	void add(int setA, int setB);

	const diff_t & differences(int setA) const ;
	const diff_t & differences(int setA,int setB) const;
	const auto & differences()  const {return _differences;} 

	const distances_t & distances(int setA) const;
	const distances_t & distances(int setA,int setB) const;
private:
	std::vector<diff_t> _differences; // vectorial distances
	std::vector<distances_t> _distances; // scalar distances
	std::map<std::pair<int,int> , int> indices2b;
	std::unordered_map<int,int> indices1b;
	geometry_t *  geo;
};

#endif