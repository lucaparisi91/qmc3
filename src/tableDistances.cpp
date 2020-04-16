#include "tableDistances.h"
#include "geometry.h"


void tableDistances::add(int setA)
{
	auto index = _differences.size();
	_differences.resize(_differences.size() + 1);
	_distances.resize(_distances.size() + 1);

	indices1b[setA] = index;
};

void tableDistances::add(int setA,int setB)
{
	auto index = _differences.size();
	_differences.resize(_differences.size() + 1);
	_distances.resize(_distances.size() + 1);

	indices2b[ std::make_pair(setA,setB)] = index;
};

const tableDistances::diff_t & tableDistances::differences(int setA) const 
{
	int index= indices1b.at(setA);
	return _differences[index];
};

const tableDistances::distances_t & tableDistances::distances(int setA) const 
{
	int index= indices1b.at(setA);
	return _distances[index];
};

const tableDistances::diff_t & tableDistances::differences(int setA, int setB) const 
{
	int index= indices2b.at(std::make_pair(setA,setB));
	return _differences[index];
};

const tableDistances::distances_t & tableDistances::distances(int setA,int setB) const 
{
	int index= indices2b.at(std::make_pair(setA,setB));

	return _distances[index];
};

void tableDistances::update(const tableDistances::states_t & states)
{
	
	for ( const auto  element : indices1b ) // updates one body distances
	{
		_differences[element.second]=geo->differences(states[element.first],{0.,0.,0.});
		_distances[element.second]=norm(_differences[element.second]);

	}

	for ( const auto  element : indices2b ) // updates one body distances
	{
		int setA=element.first.first;
		int setB= element.first.second;
		if (setA == setB)
			{
				_differences[element.second]=geo->differences(states[setA]);
			}
			else
			{
				_differences[element.second]=geo->differences(states[setA],states[setB]);

			}
		_distances[element.second]=norm(_differences[element.second]);

	}	
	
}
