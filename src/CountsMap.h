#ifndef CountsMap_h
#define CountsMap_h

#include "DataCube.h"

/** 
 * @class CountsMap
 *
 * @brief Container for binning LAT photon events over all possible
 * data dimensions including position, energy, and LAT internal 
 * degrees-of-freedom.
 *
 * @author J. Chiang
 *
 * $Header$
 */

class CountsMap : DataCube {

public:
    
    CountsMap();
    ~CountsMap();
    
private:

   //!  projection identifier, e.g., "tangent", "plat-carre", "None", etc.
   std::string m_projtype;
        
};

#endif // CountsMap_h
