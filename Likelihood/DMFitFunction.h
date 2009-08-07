/** 
 * @file DMFitFunction.h
 * @brief Declaration for the DMFit Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/DMFitFunction.h,v 1.2 2008/10/12 21:25:05 cohen Exp $
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/DMFitFunction.h,v 1.2 2008/10/12 21:25:05 cohen Exp $
 */
    
class DMFitFunction : public optimizers::Function {

public:

  DMFitFunction() : m_filename(""){
      init(1., 1., 100., 1.0, 1, 1);
   }

   /// @param norm Normalization of the function
   /// @param mass The mass of the Dark Matter particle
   /// @param bratio The branching ratio between the 2 allowed final states
   /// @param channel0 : index of the first final state
   /// @param channel1 : index of the second final state   
  DMFitFunction(double norm, double sigmav, double mass, double bratio, 
                 int channel0, int channel1) : m_filename(""){
      init(norm, sigmav, mass, bratio, channel0, channel1);
   }

   double value(optimizers::Arg&) const;

   double derivByParam(optimizers::Arg &x, const std::string &paramName) const;

   virtual Function *clone() const {
      return new DMFitFunction(*this);
   }

   void readFunction(const std::string & filename);

   const std::string & filename() const {
      return m_filename;
   }

protected:

   double integral(optimizers::Arg &, optimizers::Arg &) const {
      return 0;
   }

private:

   void init(double norm, double sigmav, double mass, double bratio, 
                 int channel0, int channel1);

   std::string m_filename;
};

} // namespace Likelihood

#endif // Likelihood_DMFitFunction_h
