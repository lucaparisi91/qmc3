#include "geometry.h"


real_t norm(real_t x,real_t y, real_t z)
{
	return std::sqrt(x*x + y*y + z*z);
}

geometryPBC::diff_t geometryPBC::differencesOneBody(const geometryPBC::particles_t & particleData, const std::array<real_t,3> & x) const
	{
		int N = particleData.dimensions()[0];
		int D=particleData.dimensions()[1];

		diff_t diffs(N,D);

		for (int i=0;i<N;i++)
		{
			for(int d=0;d<D;d++)
				{
				diffs(i,d)=difference( particleData(i,d) - x[d] ,d);
				}
		}

		return diffs;
	}


Eigen::Tensor<real_t,1> norm( const Eigen::Tensor<real_t, 2> & diffs)
{
	const int N = diffs.dimensions()[0];
	const int D= diffs.dimensions()[1];

	Eigen::Tensor<real_t,1> norms(N);
	for(int i=0;i<N;i++)
	{
		real_t tmp=0;
		for(int d=0;d<D;d++)
		{
			tmp+= diffs(i,d)*diffs(i,d);	
		}
		norms(i)=std::sqrt(tmp);
	}

	return norms;
}


geometryPBC::diff_t geometryPBC::differencesTwoBody(const geometryPBC::particles_t & particleData) const 
	{
		int N = particleData.dimensions()[0];
		int D=particleData.dimensions()[1];

		diff_t diffs( (N*(N-1))/2,D);

		int k=0;
 		for(int i=0;i<N;i++)
 		{
 			for(int j=0;j<i;j++)
 			{
 				for (int d=0;d<D;d++)
 					diffs(k,d)=difference( particleData(i,d) - particleData(j,d) , d);
 				k++;
 			}	
		}
		return diffs;
	}


 geometryPBC::diff_t geometryPBC::differencesTwoBody(const geometryPBC::particles_t & data1,const geometryPBC::particles_t & data2) const 
	{
		int N1 = data1.dimensions()[0];
		int N2 = data2.dimensions()[0];

		int D=data1.dimensions()[1];

		assert( data1.dimensions()[1] == data2.dimensions()[1]);
		diff_t diffs(N1*N2,D);

		int k=0;
 		for(int i=0;i<N1;i++)
 		{
 			for(int j=0;j<N2;j++)
 			{
 				for (int d=0;d<D;d++)
 					diffs(k,d)=difference( data1(i,d) - data2(j,d) , d);
 				k++;
 			}	
		}
		return diffs;
	}



geometry::diff_t geometry::differencesOneBody(const geometry::particles_t & data, std::initializer_list<real_t> l) const
	{
		std::array<real_t, 3> arr;
		int d=0;
		for (real_t x : l)
		 {
		 	arr[d]=x;d++;
		 }
		return differencesOneBody(data,arr);

	 }