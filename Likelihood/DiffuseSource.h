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
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools/Likelihood/src/DiffuseSource.h,v 1.1.1.1 2003/01/30 23:23:03 burnett Exp $
 */

class DiffuseSource : public Source {
    
public:
    
    DiffuseSource();
    ~DiffuseSource();
    
private:
        
};

#endif // DiffuseSource_h
