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
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools/Likelihood/src/Response.h,v 1.1.1.1 2003/01/30 23:23:03 burnett Exp $
 */

class Response : DataCube, Function, FitsNtuple {
    
public:
    
    Response();
    ~Response();
    
private:
        
};

#endif // Response_h
