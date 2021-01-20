#include <string>

class hdf5IO
{
    hdf5IO(std::string filename);

    public:

    void write( const double * data , const std::string & label, const int * dimensions ); 
    void read( const double * data , const std::string & label, const int * dimensions );

    private:


};