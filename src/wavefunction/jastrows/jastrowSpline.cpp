#include "jastrowSpline.h"


Eigen::ArrayXd getEigenArray(const json_t & j)
{
  auto vec = j.get<std::vector<real_t> >();

  Eigen::ArrayXd eigenVec;
  
  for (auto el : vec)
    {
      eigenVec << el;
    }

  return eigenVec;
  
}


jastrowSpline::jastrowSpline(const json_t & j)
  : jastrowSpline(  getEigenArray(j["coefficients"]) ,    j["step_size"].get<real_t>() , j["longDistanceConsta"].get<real_t>() ) {}



jastrowSpline::jastrowSpline(Eigen::ArrayXd coefficients_,real_t stepSize_,real_t longDistanceConstant) :
  coefficients(coefficients_),
  stepSize(stepSize_),
  inverseStepSize(1./stepSize),
  X(4),
  AX(4),
  AXd1(4),
  AXd2(4),
  A (
     (Eigen::Matrix4d() <<
      
      1/6. , 2./3 , 1./6 , 0 ,
      
      -0.5 , 0 , 0.5 , 0 ,

      0.5, -1 , 0.5 , 0 ,

      -1./6 , 0.5 , -0.5 , 1./6
      )
     .finished().transpose()
     ),
  Ad1 (
     (Eigen::Matrix4d() <<
      -0.5 , 1 , -0.5 , 0 ,
      0 , -2 , 3./2, 0,
      0.5 , 1 , -3./2 , 0,
      0, 0 , 0.5 , 0
      )
     .finished()
       ),
  Ad2 (
       (Eigen::Matrix4d() <<
	
	1 , -1 , 0 , 0 ,
	-2 , 3. , 0, 0,
	1 , -3 , 0, 0,
	0, 1, 0 , 0

	)
       .finished()
       )
{
  
  cutOff=(coefficients.size()-3  ) * stepSize;
  X.setConstant(1);
  AX.setConstant(0);
  AXd1.setConstant(0);
  AXd2.setConstant(0);

}



void jastrowSpline::evaluateDerivatives(real_t x, real_t & d0, real_t & d1, real_t & d2) 
  {
    
    if (x < cutOff)
      {
	x=x*inverseStepSize;
	int i =std::floor(x);
	x-=i;
    
	fillX(x);
    
	AX=A*X;
	AXd1=Ad1*X;
	AXd2=Ad2*X;

	
	d0=coefficients[i]*AX[0] + coefficients[i+1]*AX[1] + coefficients[i+2]*AX[2] + coefficients[i+3]*AX[3] ;
	d1=(coefficients[i]*AXd1[0] + coefficients[i+1]*AXd1[1] + coefficients[i+2]*AXd1[2] + coefficients[i+3]*AXd1[3])*inverseStepSize ;
	d2=(coefficients[i]*AXd2[0] + coefficients[i+1]*AXd2[1] + coefficients[i+2]*AXd2[2] + coefficients[i+3]*AXd2[3])*inverseStepSize*inverseStepSize ;
	
      }
    else
      {
	d0=longDistanceConstant;
	d1=0;
	d2=0;
      }

        
   
  }
