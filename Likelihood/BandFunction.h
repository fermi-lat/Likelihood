/** 
 * @file BandFunction.h
 * @brief Declaration for the Band Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/BandFunction.h,v 1.3 2015/03/03 18:05:36 jchiang Exp $
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
 */
    
class BandFunction : public optimizers::Function {

public:
 
   /// @param norm Normalization of the function
   /// @param alpha The low energy photon index
   /// @param beta The high energy photon index
   /// @param Ep The energy of the nuFnu peak (MeV)
   /// @param Scale Energy scale for power-law components (MeV)
   BandFunction(double norm=1., double alpha=-1., double beta=-2.,
                double Ep=0.1, double scale=0.1);

   virtual Function * clone() const {
      return new BandFunction(*this);
   }

protected:

   double value(const optimizers::Arg &) const;

   double derivByParamImp(const optimizers::Arg & x,
                          const std::string & paramName) const;

};

} // namespace Likelihood

#endif // Likelihood_BandFunction_h
