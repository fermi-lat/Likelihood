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
 * $Header$
 */

class FitsNtuple : public Ntuple {
    
public:
    
   FitsNtuple();
   ~FitsNtuple();

   void read_FITS_table(std::string file);
    
private:
        
};

#endif // FitsNtuple_h
