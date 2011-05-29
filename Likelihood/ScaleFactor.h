/**
 * @brief Decorator class for multiplying underlying spectral models
 * by an overall scale factor (on top of the usual normalization
 * parameter.)  Requested by DM group.
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/ScaleFactor.h,v 1.5 2009/06/19 06:36:14 jchiang Exp $
 */

#include "optimizers/Function.h"

namespace Likelihood {

class ScaleFactor : public optimizers::Function {
   
public:

   ScaleFactor();

   ScaleFactor(const optimizers::Function & spectrum,
               double scale_factor=1);

   ScaleFactor(const ScaleFactor & other);

   ScaleFactor & operator=(const ScaleFactor & rhs);

   ~ScaleFactor() throw();

   virtual double value(optimizers::Arg & x) const;
   virtual double derivByParam(optimizers::Arg & x,
                               const std::string & paramName) const;
   virtual Function * clone() const {
      return new ScaleFactor(*this);
   }

   /// Set a Parameter using a Parameter object.  This version
   /// preserves the references to the m_spectrum parameters.
   virtual void setParam(const optimizers::Parameter & param);

private:
   
   optimizers::Function * m_spectrum;

   void init(double scale_factor=1);

   void setParRefs();

};

} // namespace Likelihood
