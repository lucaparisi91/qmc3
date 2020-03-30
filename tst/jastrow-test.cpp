#include "gtest/gtest.h"
#include "parameters.h"
#include "wavefunction/jastrows/jastrow.h"
#include <vector>

TEST(jastrowTest,gradient_parameters)
{
  auto alpha=parameter("alpha",0,1);

  auto j = gaussianJastrow(1.);
  std::vector<real_t> parameters {0};
  j.registerParameter(0,alpha);
  j.addGradientParameters(1., parameters);
  
  EXPECT_EQ(parameters[0],-1.);
}