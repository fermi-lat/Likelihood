#include "../Likelihood/Function.h"

namespace Likelihood {

/** 
 * @class MyFun
 *
 * @brief a simple test function that inherits from Function
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools/Likelihood/src/MyFun.h,v 1.1 2003/02/19 01:34:34 jchiang Exp $
 */
    
class MyFun : public Function {
public:

   MyFun();
   virtual double value(double) const;
   virtual double operator()(double x) const {return value(x);};
   virtual double derivByParam(double, const std::string &paramName) const;

private:
};

} // namespace Likelihood

