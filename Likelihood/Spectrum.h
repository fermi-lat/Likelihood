#ifndef Spectrum_h
#define Spectrum_h

#include "Function.h"

/** 
 * @class Spectrum
 *
 * @brief Additive spectral components such as a power-law continuum
 * or a Gaussian emission line.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools/Likelihood/src/Spectrum.h,v 1.1.1.1 2003/01/30 23:23:03 burnett Exp $
 */

class Spectrum : public Function {
    
protected:

public:
    
   Spectrum();
   ~Spectrum();

private:
        
};

#endif // Spectrum_h
