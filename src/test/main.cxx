// test program for Likelihood

// include everything for the compiler to test
#include "Source.h"
#include "DiffuseParametric.h"
#include "DiffuseSource.h"
#include "PointSource.h"
#include "Response.h"
#include "SpectralFilter.h"
#include "Spectrum.h"
#include "Statistic.h"


int main(int argc, char * argv[]){

// make a simple test function guy that inherits from Function
    
    class MyFun: public Function {
    public:
        double value() const { return 99;}
    };

    MyFun f; 

    std::cout << "Simple test of a function derived from Function: f() = " << f() << std::endl;

    return 0;
}