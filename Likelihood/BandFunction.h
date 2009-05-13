/** 
 * @file BandFunction.h
 * @brief Declaration for the Band Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/BandFunction.h,v 1.1 2005/01/26 06:53:36 jchiang Exp $
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/BandFunction.h,v 1.1 2005/01/26 06:53:36 jchiang Exp $
 */
    
class BandFunction : public optimizers::Function {

public:

   BandFunction() {
      init(1., -1., -2., 0.1, 0.1);
   }

   /// @param norm Normalization of the function
   /// @param alpha The low energy photon index
   /// @param beta The high energy photon index
   /// @param Ep The energy of the nuFnu peak (MeV)
   /// @param Scale Energy scale for power-law components (MeV)
   BandFunction(double norm, double alpha, double beta, double Ep,
                double scale=0.1) {
      init(norm, alpha, beta, Ep, scale);
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

   void init(double norm, double alpha, double beta, double Ep,
             double scale);

};

} // namespace Likelihood

#endif // Likelihood_BandFunction_h
