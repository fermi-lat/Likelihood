#ifndef FitsNtuple_h
#define FitsNtuple_h

#include "Ntuple.h"

/** 
 * @class FitsNtuple
 *
 * @brief NTuple for FITS table data
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools/Likelihood/src/FitsNtuple.h,v 1.1.1.1 2003/01/30 23:23:03 burnett Exp $
 */

class FitsNtuple : public Ntuple {
    
public:
    
   FitsNtuple();
   ~FitsNtuple();

   void read_FITS_table(std::string file);
    
private:
        
};

#endif // FitsNtuple_h
