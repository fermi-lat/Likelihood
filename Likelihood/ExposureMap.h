#ifndef ExposureMap_h
#define ExposureMap_h

#include "DataCube.h"

/** 
 * @class ExposureMap
 *
 * @brief The integrated exposure as a function true photon direction and
 * energy.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools/Likelihood/src/ExposureMap.h,v 1.1.1.1 2003/01/30 23:23:03 burnett Exp $
 */

class ExposureMap : DataCube {
    
public:
    
    ExposureMap();
    ~ExposureMap();
    
private:

   //! projection identifier, e.g., "tangent", "plat-carre", "None", etc.
   std::string m_projtype;
                             
};

#endif // ExposureMap_h
