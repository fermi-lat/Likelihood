/** @file MyFun.h
 * @brief Test function declaration.
 * @author J. Chiang
 *
 * $Header$
 */

#include "Likelihood/Function.h"
#include "Likelihood/Arg.h"

namespace Likelihood {

/** 
 * @class MyFun
 *
 * @brief A simple test function that inherits from Function
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/MyFun.h,v 1.5 2003/03/16 21:53:26 jchiang Exp $
 */
    
class MyFun : public Function {
public:

   MyFun();
   ~MyFun(){}

   double value(Arg &) const;

   double derivByParam(Arg &x, const std::string &paramName) const;

private:

};

} // namespace Likelihood

