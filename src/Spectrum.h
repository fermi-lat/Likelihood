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
 * $Header$
 */

class Spectrum : public Function {
    
protected:

public:
    
   Spectrum();
   ~Spectrum();

private:
        
};

#endif // Spectrum_h
