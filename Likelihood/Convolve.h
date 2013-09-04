/**
 * @file Convolve.h
 * @brief Static methods to perform 1D and 2D convolutions using the
 *        FFTW library
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/Convolve.h,v 1.1 2005/05/23 05:51:25 jchiang Exp $
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

   static std::vector<float> 
   convolve(const std::vector<float> & signal,
            const std::vector<float> & psf);

   static std::vector< std::vector<float> > 
   convolve2d(const std::vector< std::vector<float> > & signal,
              const std::vector< std::vector<float> > & psf);

};

} // namespace Likelihood

#endif // Likelihood_Convolve_h
