/** 
 * @file DMFitFunction.h
 * @brief Declaration for the DMFit Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/DMFitFunction.h,v 1.4 2011/11/29 16:41:17 cohen Exp $
 */

#ifndef Likelihood_DMFitFunction_h
#define Likelihood_DMFitFunction_h

#include "optimizers/Function.h"
#include "optimizers/Arg.h"

#include <string>

namespace Likelihood {

/** 
 * @class DMFitFunction
 *
 * @brief A DMFit function for modeling Dark Matter spectra.
 *
 * @author J. Cohen-Tanugi, based on the DMFit package by S. Profumo and T. Jeltema
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/DMFitFunction.h,v 1.4 2011/11/29 16:41:17 cohen Exp $
 */
    
class DMFitFunction : public optimizers::Function {

public:

   /// @param norm Normalization of the function
   /// @param mass The mass of the Dark Matter particle
   /// @param bratio The branching ratio between the 2 allowed final states
   /// @param channel0 : index of the first final state
   /// @param channel1 : index of the second final state   
   DMFitFunction(double norm=1., double sigmav=1., double mass=100., 
                 double bratio=1., int channel0=1, int channel1=1);

   virtual Function *clone() const {
      return new DMFitFunction(*this);
   }

   void readFunction(const std::string & filename);

   const std::string & filename() const {
      return m_filename;
   }

protected:

   virtual double value(optimizers::Arg&) const;

   virtual double derivByParamImp(optimizers::Arg & x, 
                                  const std::string & paramName) const;

   virtual double integral(optimizers::Arg &, optimizers::Arg &) const {
      return 0;
   }

private:

   std::string m_filename;

   double m_8pi;

};

} // namespace Likelihood

#endif // Likelihood_DMFitFunction_h
