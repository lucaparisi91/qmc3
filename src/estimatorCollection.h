#ifndef ESTIMATORS_COLLECTION_H
#define ESTIMATORS_COLLECTION_H


#include <vector>
#include <memory>

class estimatorBase;

class estimatorCollection
{
public:
	estimatorCollection(){};
	auto begin() {return _estimators.begin();}
	const auto cbegin() const {return _estimators.begin();}

	auto end() {return _estimators.end();}
	const auto cend() const {return _estimators.end();}

        void push_back(estimatorBase* est){_estimators.emplace_back(est);}

	auto & operator[](size_t i) {return _estimators[i];}
	virtual void clear();
	virtual void dump();
        virtual void accumulateMPI(int root=0);
  
private:
  std::vector<std::unique_ptr<estimatorBase> > _estimators;
};

#endif
