#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "traits.h"
#include <tuple>
#include <cmath>
#include <initializer_list>
#include <unsupported/Eigen/CXX11/Tensor>



class geometry
{
public:
	using particles_t = ::state_t;
	using diff_t =  ::difference_t;


	virtual diff_t differencesTwoBody(const particles_t & data) const = 0; // distances between undistinguishible particles

	virtual diff_t differencesTwoBody(const particles_t & data1, const particles_t & data2) const = 0; // distances between distinguishiblle particles


	virtual diff_t differencesOneBody(const particles_t & data , const std::array<real_t,3> & x) const = 0;

	virtual diff_t differencesOneBody(const particles_t & data, std::initializer_list<real_t> l) const;

	virtual diff_t differences(const particles_t & data, std::initializer_list<real_t> l) const { return differencesOneBody(data,l);}

	virtual diff_t differences(const particles_t & data) const { return differencesTwoBody(data);}

	virtual diff_t differences(const particles_t & data1, const particles_t & data2 ) const { return differencesTwoBody(data1,data2);}

};


real_t norm(real_t x,real_t y, real_t z);

Eigen::Tensor<real_t,1> norm( const Eigen::Tensor<real_t, 2> & diffs);


class geometryOBC
{
public:
   void difference(real_t diffx,real_t diffy,real_t diffz,real_t & diffNewx, real_t & diffNewy,real_t & diffNewz) {diffNewx=diffx;diffNewy=diffy;diffNewz=diffz;}
   std::tuple<double,double,double> difference(real_t diffx,real_t diffy, real_t diffz) { return std::make_tuple(diffx,diffy,diffz);}
private:
};



class geometryPBC : public geometry
{
public:
	geometryPBC(real_t lBoxx_,real_t lBoxy_,real_t lBoxz_) : lBox{ lBoxx_,lBoxy_,lBoxz_},lBoxInverse{1./lBoxx_,1./lBoxy_,1./lBoxz_} {} ;

	real_t difference(real_t t, int i ) const {return ( t - std::round(t*lBoxInverse[i] )*lBox[i]);}

	auto difference(real_t diffx,real_t diffy,real_t diffz) const { return std::make_tuple( difference(diffx,0) , difference(diffy,1) , difference(diffz,1)   );}

	virtual diff_t differencesOneBody(const particles_t & particleData, const std::array<real_t,3> & x) const;

	virtual diff_t differencesTwoBody(const particles_t & particleData) const ;

	virtual diff_t differencesTwoBody(const particles_t & data1,const particles_t & data2) const ;

private:
	real_t lBox [3];
	real_t lBoxInverse [3];

};

#endif