#ifndef Response_h
#define Response_h

#include "DataCube.h"
#include "Function.h"
#include "FitsNtuple.h"

/** 
 * @class Response
 *
 * @brief The LAT instrumental response.  Multiple representations
 * may be required:
 *
 * DataCube: in case the standard decomposition into effective area, energy
 * dispersion and point-spread function does not apply
 *
 * Function: for use with log-likelihood sums and integrals
 *
 * FitsNtuple: for the straw-man CALDB representation
 *
 * @author J. Chiang
 *    
 * $Header$
 */

class Response : DataCube, Function, FitsNtuple {
    
public:
    
    Response();
    ~Response();
    
private:
        
};

#endif // Response_h
