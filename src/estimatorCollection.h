#include <vector>

class estimatorBase;

class estimatorCollection
{
public:
	estimatorCollection(){};
	auto begin() {return _estimators.begin();}
	const auto cbegin() const {return _estimators.begin();}

	auto end() {return _estimators.end();}
	const auto cend() const {return _estimators.end();}

	void push_back(estimatorBase* est){_estimators.push_back(est);}

	auto & operator[](size_t i) {return _estimators[i];}
	virtual void clear();
	virtual void dump();
	
private:
	std::vector<estimatorBase*> _estimators;
};
