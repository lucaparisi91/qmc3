#include "potential.h"

class harmonicPotential
{
public:
	harmonicPotential(real_t omegax,real_t omegay,real_t omegaz);

	real_t operator()(const states_t & states, const differences_t & differences, const distances_t & distances ) const;

private:
	std::array<real_t, 3> frequencies;
}