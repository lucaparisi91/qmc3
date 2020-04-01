#include "traits.h"
#include "qmcExceptions.h"

class tableDistances;

class potential
{
	/*
	Represents the logarithm of the wavefunction. Supports wavefunctions acting on at most two different sets
	*/
public:
	using states_t=::states_t;
	using state_t = ::state_t;

	potential(const geometry_t & geo_ );

	virtual real_t operator()(const states_t & state ) = 0;
	virtual real_t operator()(const states_t & state , const tableDistances & tab) = 0;

	const auto & getGeometry() {return *geo;}

private:
	const geometry_t * geo;
};

class potential1b : public potential
{
public:

	potential1b(const geometry_t & geo, int setA_); 

	virtual real_t operator()(const states_t & states ){return (*this)(states[_setA]);}
	virtual real_t operator()(const states_t & state,const tableDistances & tab ) ;



	virtual real_t operator()(const state_t & state);

	virtual real_t operator()(const state_t & state,const distance_t & dis){throw missingImplementation("Evaluation on a single state from distances in one body potential.");}

	const int & setA() const {return _setA;}

private:
	int _setA;
	distance_t distances;
	difference_t differences;
};

class harmonicPotential :  public potential1b
{
public:
	using potential1b::operator();

	harmonicPotential(const geometry_t & geo,real_t freq , int setA ) ;

	virtual real_t operator()(const state_t & state,const distance_t & dis) override;
private:

	real_t omega;
};