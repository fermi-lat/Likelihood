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
 * $Header$
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
