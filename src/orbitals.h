#include "traits.h"

class sinOrbital
{
public:
  using value_t = real_t;
  
  sinOrbital(size_t nx,size_t ny,size_t nz,real_t lBox_);
  
  real_t operator()(real_t x,real_t y,real_t z) const { return std::sin(k[0]*x + k[1]*y + k[2]*z + delta) ; }

  real_t firstDerivative(real_t x1,real_t x2,real_t x3,size_t id) const {return k[id]*std::cos(k[0]*x1 + k[1]*x2 + k[2]*x3 + delta)  ;}
  
  real_t secondDerivative(real_t x1,real_t x2,real_t x3,size_t id) const {return -k[id]*k[id]*std::sin(k[0]*x1 + k[1]*x2 + k[2]*x3 + delta)  ;}
  
  real_t energy() const {return (k[0]*k[0] + k[1]*k[1] + k[2]*k[2] )/2.;}
  
private:
  real_t delta;
  real_t lBox;
  std::vector<int> ns;
  std::vector<int> k;
};


template<class orbital_t>
class orbitalSet
{
private:
  using matrix_t = Tensor::Eigen<value_t,2> ;
  
  void evaluateMatrix(const state_t & state,matrix_t & matrix) ; // store in a matrix all the determinants
  
public:
  std::vector<orbital_t> orbitals;
};
