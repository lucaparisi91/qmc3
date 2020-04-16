#ifndef SLATERS_H
#define SLATERS_H


#include "traits.h"
#include <map>
#include <unordered_map>

class orbitalSetBase;

double getLog(double a,int & sign);
std::complex<double> getLog(std::complex<double> a,int & sign);

class tableSlaters
{
public:
  using matrix_t = Eigen::MatrixXd;
  using states_t = ::states_t;
  
  tableSlaters(){};
  
  void update(const states_t & states); // update all registered distances
  
  void add(int setA ,orbitalSetBase * orbitalsSet); 

  //void add(int setA, int setB,orbitalSetBase * orbitalSet);
  
  const matrix_t & slaterMatrix(int setA) const { return slaterMatrices[indices1b.at(setA) ]; }
  const matrix_t & slaterMatrixInverse(int setA) const { return slaterMatricesInverse[indices1b.at(setA) ]; }
  
  
  //const slaterMatrix & slaterMatrix(int setA,int setB) const;
  
  real_t logDeterminant(int setA) const {return logDeterminants[indices1b.at(setA)];}
  size_t size() {return orbitalSets.size();}
private:
  const geometry_t *  geo;
  std::vector<orbitalSetBase* > orbitalSets;
  std::vector<matrix_t> slaterMatrices;
  std::vector<matrix_t> slaterMatricesInverse;
  std::unordered_map<int,int> indices1b;
  std::vector<real_t> logDeterminants;
  std::vector<int> signs;
  Eigen::PartialPivLU<matrix_t> lud;  
  //std::map<std::pair<int,int> , int> indices2b;
};

#endif
