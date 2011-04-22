/** 
 * @file LogParabola.h
 * @brief Declaration for the LogParabola class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/LogParabola.h,v 1.1 2005/07/18 22:54:58 jchiang Exp $
 */

#ifndef Likelihood_LogParabola_h
#define Likelihood_LogParabola_h

#include "optimizers/Function.h"
#include "optimizers/Arg.h"

namespace Likelihood {

/** 
 * @class LogParabola
 *
 * @brief A log-parabolic function for modeling Blazar SED components.
 *
 */
    
class LogParabola : public optimizers::Function {

public:

   LogParabola() {
      init(1., 1., 2., 100);
   }

   /// @brief The log-parabolic form is given by
   /// \f$K (E/E_1)^{-(\alpha + \beta*\log(E/E_1))}\f$
   /// @param norm \f$k\f$ in the formula
   /// @param alpha 
   /// @param beta 
   /// @param Eb \f$E_1\f$ in the formula
   LogParabola(double norm, double alpha, double beta, double Eb) {
      init(norm, alpha, beta, Eb);
   }

   double value(optimizers::Arg&) const;

   double derivByParam(optimizers::Arg &x, const std::string &paramName) const;

   virtual Function *clone() const {
      return new LogParabola(*this);
   }

protected:

   double integral(optimizers::Arg &, optimizers::Arg &) const {
      return 0;
   }

private:

   void init(double norm, double alpha, double beta, double Eb);

};

} // namespace Likelihood

#endif // Likelihood_LogParabola_h
