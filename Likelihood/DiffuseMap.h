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
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools/Likelihood/src/DiffuseMap.h,v 1.1.1.1 2003/01/30 23:23:02 burnett Exp $
 */

class DiffuseMap : DiffuseSource, DataCube {
    
public:
    
    DiffuseMap();
    ~DiffuseMap();
    
private:

};

#endif // DiffuseMap_h
