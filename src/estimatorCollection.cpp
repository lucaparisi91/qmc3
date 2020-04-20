#include "estimatorCollection.h"
#include "estimators.h"

void estimatorCollection::clear()
	{
		for (auto est : _estimators)
		{
			est->clear();
		}
	}


void estimatorCollection::dump()
		{
			for (auto est : _estimators)
			{
				est->dump();
			}
		}


void estimatorCollection::accumulateMPI(int root)
		{
			for (auto est : _estimators)
			{
			  est->accumulateMPI(root);
			}
		}
