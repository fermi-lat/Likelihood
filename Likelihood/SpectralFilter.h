#ifndef SpectralFilter_h
#define SpectralFilter_h

#include "Function.h"

/** 
 * @class SpectralFilter
 *
 * @brief Multiplicative spectral components such as intergalactic 
 * infrared absorption or discrete absorption features due to exotic 
 * matter along the line-of-sight.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools/Likelihood/src/SpectralFilter.h,v 1.1.1.1 2003/01/30 23:23:03 burnett Exp $
 */

class SpectralFilter : Function {
    
public:
    
   SpectralFilter();
   ~SpectralFilter();

private:
        
};

#endif // SpectralFilter_h
