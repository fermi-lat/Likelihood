/** 
 * @file dArg.h
 * @brief Declaration of dArg class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/dArg.h,v 1.4 2003/03/17 00:53:43 jchiang Exp $
 */

#ifndef dArg_h
#define dArg_h

#include "Likelihood/Arg.h"

namespace Likelihood {

/** 
 * @class dArg
 *
 * @brief Concrete Arg subclass for encapsulating data of type double.
 *
 * @authors J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/dArg.h,v 1.4 2003/03/17 00:53:43 jchiang Exp $
 */

class dArg : public Arg{
    
public:
   
   dArg(double x) : m_val(x) {}
   virtual ~dArg() {}

   double getValue() {return m_val;}

private:

   double m_val;

};

} // namespace Likelihood

#endif // dArg_h
