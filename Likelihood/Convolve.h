/**
 * @file Convolve.h
 * @brief Static methods to perform 1D and 2D convolutions using the
 *        FFTW library
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef Likelihood_Convolve_h
#define Likelihood_Convolve_h

#include <vector>

namespace Likelihood {

class Convolve {

public:

   static std::vector<double> 
   convolve(const std::vector<double> & signal,
            const std::vector<double> & psf);

   static std::vector< std::vector<double> > 
   convolve2d(const std::vector< std::vector<double> > & signal,
              const std::vector< std::vector<double> > & psf);

};

} // namespace Likelihood

#endif // Likelihood_Convolve_h
