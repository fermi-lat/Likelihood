/** 
 * @file dArg.h
 * @brief Declaration of dArg class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/dArg.h,v 1.5 2003/05/29 00:29:39 jchiang Exp $
 */

#ifndef Likelihood_dArg_h
#define Likelihood_dArg_h

#include "Likelihood/Arg.h"

namespace Likelihood {

/** 
 * @class dArg
 *
 * @brief Concrete Arg subclass for encapsulating data of type double.
 *
 * @authors J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/dArg.h,v 1.5 2003/05/29 00:29:39 jchiang Exp $
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

#endif // Likelihood_dArg_h
