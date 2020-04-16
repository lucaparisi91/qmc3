#include "traits.h"

class orbitalSetBase
{
public:
  using matrix_t = Eigen::MatrixXd;
  
  virtual void storeEvaluate(const state_t & state,matrix_t & mat) const=0;
  virtual real_t energy() const =0;
};


template<class orbital_t>
class orbitalSet : public orbitalSetBase
{
public:
  virtual void storeEvaluate(const state_t & state,matrix_t & matrix) const; // store in a matrix all orbitals

  auto & getOrbitals() {return orbitals;}
  const auto & getOrbitals() const{return orbitals;}
  
  real_t energy() const 
  {
    real_t e=0;
    for ( auto & orbital : orbitals)
      {
	e+=orbital.energy();
      }
    return e;
  }
private:
  std::vector<orbital_t> orbitals;
};

#if DIMENSIONS == 3




class sinOrbital
{
public:
  using value_t = real_t;
  
  sinOrbital(int nx,int ny,int nz,real_t lBox_);
  
  real_t operator()(real_t x,real_t y,real_t z) const { return std::sin(k[0]*x + k[1]*y + k[2]*z + delta) ; }
  
  inline void evaluateDerivatives(real_t x,real_t y, real_t z, real_t & dx, real_t & dy, real_t & dz ,real_t & laplacian ) const
  {
    real_t tmp = std::cos(k[0]*x + k[1]*y +k[2]*z + delta );
    real_t tmp1 = std::sin(k[0]*x + k[1]*y + k[2]*z + delta);
    dx=k[0]*tmp;
    dy=k[1]*tmp;
    dz=k[2]*tmp;
    laplacian=-(k[0]*k[0]+ k[1]*k[1]+ k[2]*k[2])*tmp1;
  }

  virtual real_t energy() const  {return 0.5*(k[0]*k[0] + k[1]*k[1] + k[2]*k[2] );}
  
private:
  real_t delta;
  real_t lBox;
  std::vector<int> ns;
  std::vector<real_t> k;
  
};

#include "wavefunction/shell.h"
template<class orbital_t>
void fillFermiSea(std::vector<orbital_t> & orbitals,int N,double lBox)
{
  int n=0;
  shellStructure sh(N);
  
    for(int i=0;i<sh.nShells() and n<N;i++)
    {
      for(int j=0;j<sh[i].capacity() and n < N;j++)
	{
	  auto indices=sh[i][j];
	  orbitals.push_back(orbital_t(std::get<0>(indices),std::get<1>(indices),std::get<2>(indices),lBox));
	  n++;
	}
    }    
};


#endif
