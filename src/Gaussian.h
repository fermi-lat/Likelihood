/** 
 * @file Gaussian.h
 * @brief Gaussian class declaration
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Gaussian.h,v 1.10 2003/07/19 04:38:03 jchiang Exp $
 */

#ifndef Likelihood_Gaussian_h
#define Likelihood_Gaussian_h

#include "optimizers/Function.h"
namespace optimizers {
   class Arg;
}

namespace Likelihood {
/** 
 * @class Gaussian
 *
 * @brief A 1D Gaussian function
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Gaussian.h,v 1.10 2003/07/19 04:38:03 jchiang Exp $
 */
    
class Gaussian : public optimizers::Function {
public:

   Gaussian(){init(0, 0, 1);}
   Gaussian(double Prefactor, double Mean, double Sigma)
      {init(Prefactor, Mean, Sigma);}

   double value(optimizers::Arg &) const;

   double derivByParam(optimizers::Arg &, const std::string &paramName) const
      throw(optimizers::ParameterNotFound);

   double integral(optimizers::Arg &xmin, optimizers::Arg &xmax) const;

   virtual optimizers::Function *clone() const {
      return new Gaussian(*this);
   }

private:

   void init(double Prefactor, double Mean, double Sigma);
   double erfcc(double x) const;

};

} // namespace Likelihood

#endif // Likelihood_Gaussian_h
