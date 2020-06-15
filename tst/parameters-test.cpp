#include "gtest/gtest.h"
#include "wavefunction/jastrows/jastrow.h"
#include <vector>

TEST(parametersTest,init)
{

  optimizationParameters parameters;
  bool status=parameters.addParameter("alpha",1);
  ASSERT_EQ(status,true);
  
  status=parameters.addParameter("alpha",1);
  ASSERT_EQ(status,false);
  
  optimizationParameter alpha("alpha",1,0);
  auto param = parameters.mapParameter(alpha,"alpha");
  
  
  gaussianJastrow J(1.0);
  
  
  std::vector<real_t> gradPs(1,0);
  
  J.addGradientParameter(1.,alpha,gradPs.begin(),gradPs.end() );
  
  ASSERT_EQ(gradPs[0],-1.);

  gradPs={0};
  
  J.addGradientParameter(1., param,gradPs);
  
  ASSERT_EQ(gradPs[0],-1.);
  
}
