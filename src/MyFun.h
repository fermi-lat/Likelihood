#include "Likelihood/Function.h"
#include "Likelihood/Arg.h"

namespace Likelihood {

/** 
 * @class MyFun
 *
 * @brief a simple test function that inherits from Function
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/MyFun.h,v 1.3 2003/03/04 17:45:32 jchiang Exp $
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

