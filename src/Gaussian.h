/** 
 * @file Gaussian.h
 * @brief Gaussian class declaration
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Gaussian.h,v 1.9 2003/05/29 20:10:46 jchiang Exp $
 */

#ifndef Likelihood_Gaussian_h
#define Likelihood_Gaussian_h

#include "Likelihood/Function.h"
#include "Likelihood/Arg.h"

namespace Likelihood {
/** 
 * @class Gaussian
 *
 * @brief A 1D Gaussian function
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Gaussian.h,v 1.9 2003/05/29 20:10:46 jchiang Exp $
 */
    
class Gaussian : public Function {
public:

   Gaussian(){init(0, 0, 1);}
   Gaussian(double Prefactor, double Mean, double Sigma)
      {init(Prefactor, Mean, Sigma);}

   double value(Arg &) const;

   double derivByParam(Arg &, const std::string &paramName) const
      throw(ParameterNotFound);

   double integral(Arg &xmin, Arg &xmax) const;

   virtual Function *clone() const {
      return new Gaussian(*this);
   }

private:

   void init(double Prefactor, double Mean, double Sigma);
   double erfcc(double x) const;

};

} // namespace Likelihood

#endif // Likelihood_Gaussian_h
