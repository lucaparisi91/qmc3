#include "initializer.h"
#include "wavefunction/productWavefunction.h"
#include "tableDistances.h"
#include "qmcExceptions.h"
#include "walkers.h"
#include "wavefunction/wavefunction.h"


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

void initializer::registerSlaterOrbitals(tableSlaters & tab,const productWavefunction & psi)
{
  for (const auto & wave : psi.waves() )
    {
      const auto  sets = wave->sets();
      const auto  orbitals = wave->orbitals();

      if (orbitals.size() == 1 and sets.size() == 1 )
	{
	  tab.add(sets[0],orbitals[0]);
	}
    }
}

void initializer::initialize(walker & w, const states_t & states ,  productWavefunction & psi)
{
	w.getStates()=states;

	registerDistances(w.getTableDistances(),psi);
	registerSlaterOrbitals(w.getTableSlaters(),psi);
	w.getTableDistances().update(states);
	w.getTableSlaters().update(states);
	psi.evaluateDerivatives( w);
	
}

void initializer::initialize(dmcWalker & w, const states_t & states ,  productWavefunction & psi,energy & ob)
{
  w.getStates()=states;
  registerDistances(w.getTableDistances(),psi);
  registerSlaterOrbitals(w.getTableSlaters(),psi);
  updateForceGradientEnergy(w,psi,ob);
}

void initializer::initialize(walkerContainer<dmcWalker> & ws, const std::vector<states_t> & states ,  productWavefunction & psi,energy & ob)
{
  ws.resize(states.size());
  
  for (int i=0;i<states.size();i++)
    {
      initialize(ws[i],states[i],psi,ob);
    }
}


  
