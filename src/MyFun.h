#include "../Likelihood/Function.h"
#include "../Likelihood/Arg.h"

namespace Likelihood {

/** 
 * @class MyFun
 *
 * @brief a simple test function that inherits from Function
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/MyFun.h,v 1.2 2003/02/23 22:09:48 jchiang Exp $
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

