#ifndef DiffuseMap_h
#define DiffuseMap_h

#include "DiffuseSource.h;
#include "DataCube.h"

/** 
 * @class DiffuseMap
 *
 * @brief Container for a diffuse emission model based on pre-computed
 * FITS images.  A standard example would be a Galactic diffuse emission 
 * model.
 *
 * @author J. Chiang
 *    
 * $Header$
 */

class DiffuseMap : DiffuseSource, DataCube {
    
public:
    
    DiffuseMap();
    ~DiffuseMap();
    
private:

};

#endif // DiffuseMap_h
