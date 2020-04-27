#include "driver.h"
#include "energy.h"

#include "walkers.h"

class realScalarEstimator;
class branchingControl;


class acceptRejectPolicy
{
public:
	virtual bool accept(
		mover & move,
		const dmcWalker & new_walker,
		const dmcWalker & old_walker,
		wavefunction_t & wavefunction,
			    randomGenerator_t & generator)=0;
  virtual real_t getAcceptanceRatio() const =0;

  virtual void clear()=0;

  virtual void accumulateMPI(int root)=0;
};

class noMetropolisPolicy : public acceptRejectPolicy
{
public:
  virtual bool accept(
		mover & move,
		const dmcWalker & new_walker,
		const dmcWalker & old_walker,
		wavefunction_t & wavefunction,
		randomGenerator_t & generator);
  virtual real_t getAcceptanceRatio() const ;
  virtual void clear() {};
  virtual void accumulateMPI(int root) {};
private:
  metropolis metropolisSampler;
};


class metropolisPolicy : public acceptRejectPolicy
{
public:
  virtual bool accept(
		mover & move,
		const dmcWalker & new_walker,
		const dmcWalker & old_walker,
		wavefunction_t & wavefunction,
		randomGenerator_t & generator
		      );
  virtual real_t getAcceptanceRatio() const ;
  virtual void clear() ;
  virtual void accumulateMPI(int root);

private:
  metropolis metropolisSampler;
  
};

class dmcDriver : public driver
{
	using walker_t =dmcWalker;
public:
  
  dmcDriver(wavefunction_t * wave,potential_t * pot,real_t timeStep,size_t nWalkers);
  
  virtual void run( const std::vector<states_t> &states , size_t nBlocks );

  virtual void step();

  void update( dmcWalker & wNew,const dmcWalker & wOld);
  
  virtual void out();
  virtual void accumulate();

  virtual void disableBranching();
  
private:
  walkerContainer<dmcWalker> current_walkers;
  walkerContainer<dmcWalker> old_walkers;
  std::unique_ptr<mover> dmcMover;
  energy energyOb;
  energyFromWalker energyAccFromWalker;
  std::unique_ptr<realScalarEstimator> energyEst;
  std::unique_ptr<acceptRejectPolicy> accepter;
  std::unique_ptr< branchingControl> brancher;
  std::unique_ptr<pTools::walkerDistribution> walkerLoadBalancer;
  bool performBranching;
};
