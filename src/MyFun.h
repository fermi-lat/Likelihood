/** 
 * @file MyFun.h
 * @brief Test function declaration.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/MyFun.h,v 1.6 2003/03/17 00:53:44 jchiang Exp $
 */

#include "Likelihood/Function.h"

namespace Likelihood {

class Arg;

/** 
 * @class MyFun
 *
 * @brief A simple test function that inherits from Function
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/MyFun.h,v 1.6 2003/03/17 00:53:44 jchiang Exp $
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

