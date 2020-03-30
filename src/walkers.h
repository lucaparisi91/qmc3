class dmcWalker
{
public:
	dmcWalker();
private:
	real_t e; // stores the energy of the current configuration
	states_t states; // a vector of particle data
	real_t waveValue; // log of the wavefunction
	grades_t gradients; //  vector of the gradients of the wavefunction for each sets
};

class dmcWalkerDistanceTable : public dmcWalker
{
public:
	dmcWalkerDistanceTable();
private:
	tableDistances distanceTable;
};