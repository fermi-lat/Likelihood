/** 
 * @file BandFunction.h
 * @brief Declaration for the Band Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/optimizers/src/BandFunction.h,v 1.2 2004/12/22 06:00:40 jchiang Exp $
 */

#ifndef Likelihood_BandFunction_h
#define Likelihood_BandFunction_h

#include "optimizers/Function.h"
#include "optimizers/Arg.h"

namespace Likelihood {

/** 
 * @class BandFunction
 *
 * @brief A Band function for modeling GRB spectra.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/BandFunction.h,v 1.2 2004/12/22 06:00:40 jchiang Exp $
 */
    
class BandFunction : public optimizers::Function {

public:

   BandFunction() {
      init(1., -1., -2., 0.1);
   }

   /// @param norm Normalization of the function
   /// @param alpha The low energy photon index
   /// @param beta The high energy photon index
   /// @param Ep The energy of the nuFnu peak (MeV)
   BandFunction(double norm, double alpha, double beta, double Ep) {
      init(norm, alpha, beta, Ep);
   }

   double value(optimizers::Arg&) const;

   double derivByParam(optimizers::Arg &x, const std::string &paramName) const;

   virtual Function *clone() const {
      return new BandFunction(*this);
   }

protected:

   double integral(optimizers::Arg &, optimizers::Arg &) const {
      return 0;
   }

private:

   void init(double norm, double alpha, double beta, double Ep);

};

} // namespace Likelihood

#endif // Likelihood_BandFunction_h
