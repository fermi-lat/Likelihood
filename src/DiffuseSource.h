#ifndef DiffuseSource_h
#define DiffuseSource_h

#include "Source.h"

/** 
 * @class DiffuseSource
 *
 * @brief Diffuse sources of gamma-ray emission, e.g., Galactic and
 * extragalactic diffuse emission, as well as discrete sources such
 * as the LMC, the Moon, the Sun, SNRs, etc.
 *
 * @author J. Chiang
 *    
 * $Header$
 */

class DiffuseSource : public Source {
    
public:
    
    DiffuseSource();
    ~DiffuseSource();
    
private:
        
};

#endif // DiffuseSource_h
