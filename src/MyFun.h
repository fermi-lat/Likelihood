#include "../Likelihood/Function.h"

namespace Likelihood {

/** 
 * @class MyFun
 *
 * @brief a simple test function that inherits from Function
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools/Likelihood/src/Function.h,v 1.1.1.1 2003/01/30 23:23:03 burnett Exp $
 */
    
class MyFun : public Function {
public:

   MyFun(){setMaxNumParams(3);};
   virtual double value(const double) const;
   virtual double operator()(const double x) const {return value(x);};
   virtual double derivByParam(const double, 
			       const std::string paramName) const;
   virtual std::vector<double> getDerivs(const double) const;

private:
};

} // namespace Likelihood

