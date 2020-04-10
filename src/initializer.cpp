#include "initializer.h"
#include "wavefunction/productWavefunction.h"
#include "tableDistances.h"
#include "qmcExceptions.h"
#include "walkers.h"

void initializer::registerDistances(tableDistances & tab,const wavefunction & wave)
{

	const auto & sets = wave.sets();
	tab.setGeometry(wave.getGeometry());
	

	if (sets.size() == 1)
	{
		tab.add(sets[0]);
	}
	else
	{
		throw missingImplementation("Register distances: more then 1b wavefunctions.");
	}
};

void initializer::registerDistances(tableDistances & tab,const productWavefunction & waves)
{
	
	for (const auto & wave : waves.waves() )
		registerDistances(tab,*wave);
};


void initializer::initialize(walker & w, const states_t & states ,  productWavefunction & psi)
{

	w.getStates()=states;

	registerDistances(w.getTableDistances(),psi);
	w.getTableDistances().update(states);

	real_t lap;	

	psi.evaluateDerivatives(w.getStates(), w.getGradients(), w.getLogWave(),w.getLaplacianLog(),w.getTableDistances());

}
