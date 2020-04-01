template<class T>
class scalarAccumulator
{
public:
	using value_t = T;
	scalarAccumulator(){}

	void operator+=(value_t e){sum+=e;n+=1;};

	value_t average() const {return sum/n;}

private:
	value_t sum;
	size_t n;
};

template<class T>
class vectorAccumulator
{
	using value_t = T::value_t;
public:
	vectorAccumulator(size_t size);

	void accumulate(value_t T, size_t i);

	T average() const;

	const T & sums() const {return vec;}

	size_t size(){return vec.size();}
	
private:
	T vec;
	std::vector<size_t> n;
};