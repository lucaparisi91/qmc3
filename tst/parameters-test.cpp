#include "gtest/gtest.h"
#include "parameters.h"
#include <vector>

TEST(parametersTest,init)
{
  auto alpha=parameter("alpha",0,1);
  std::vector<real_t> parameters{0.};

  *(alpha.begin(parameters))=2.;
  ASSERT_EQ(parameters[0],2.);
  
}