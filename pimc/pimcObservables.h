#include "tools.h"
#include "pimcConfigurations.h"
#include "action.h"
#include "accumulators.h"


namespace pimc
{



class scalarEstimator
{
    public:


    virtual Real operator()(configurations_t & configurations, firstOrderAction & S) = 0;
    


};

class observable
{
public:

    virtual void accumulate(configurations_t & configurations, firstOrderAction & S)=0;

    virtual void out(size_t iteration)=0;

    virtual void clear()=0;

};

class scalarObservable : public observable
{
public:

    scalarObservable(std::shared_ptr<scalarEstimator> ob_ , std::string label_) : label(label_),filename(label_ + ".dat"),delim(" "){
        ob=ob_;
        f.open(filename,std::fstream::app);
    }

    virtual void accumulate(configurations_t & configurations, firstOrderAction & S)
    {
        acc+=(*ob)(configurations,S);

    }

    virtual void out(size_t iteration)
    {
        if ( acc.getWeight() != 0 )
        {
            f << iteration << delim <<  acc.average() << std::endl;
        }
    }



    virtual void clear()
    {
        acc.clear();
    }

    ~scalarObservable()
    {
        f.close();
    }

    Real average()
    {
        return acc.average();
    }

    Real weight()
    {
        return acc.getWeight();
    }


    private:

    std::shared_ptr<scalarEstimator> ob;
    scalarAccumulator<Real> acc;
    std::ofstream f;
    std::string filename;
    std::string label;
    std::string delim;
};


class thermodynamicEnergyEstimator : public scalarEstimator 
{
    public:
    thermodynamicEnergyEstimator(){}
    virtual Real operator()(configurations_t & configurations, firstOrderAction & S);
};


class virialEnergyEstimator : public scalarEstimator
{
    public:
    virialEnergyEstimator(int nMax, int MMax) : buffer(nMax,getDimensions(),MMax),rC(nMax,getDimensions(),MMax) {}
    Real operator()(configurations_t & configurations, firstOrderAction & S);
    
    private:
    Eigen::Tensor<Real,3> buffer;
    Eigen::Tensor<Real,3> rC;
};




}