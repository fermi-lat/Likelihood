/** @file SrcArg.h
 * @brief Declaration of SrcArg class
 * $Header:
 */

#ifndef SrcArg_h
#define SrcArg_h

#include "Likelihood/Arg.h"
#include "Likelihood/Source.h"

namespace Likelihood {

/** 
 * @class SrcArg
 *
 * @brief Concrete Arg subclass for encapsulating data of type Source.
 *
 * @authors J. Chiang
 *    
 * $Header: */

class SrcArg : public Arg {
    
public:
   
   SrcArg(Source *src) : m_val(src) {}
   virtual ~SrcArg() {}

   Source *getValue() const {return m_val;}

private:

   Source *m_val;

};

} // namespace Likelihood

#endif // SrcArg_h
