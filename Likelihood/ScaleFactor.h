/**
 * @brief Decorator class for multiplying underlying spectral models
 * by an overall scale factor (on top of the usual normalization
 * parameter.)  Requested by DM group.
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/ScaleFactor.h,v 1.1 2011/05/29 17:53:06 jchiang Exp $
 */

#include "optimizers/Function.h"

namespace Likelihood {

class ScaleFactor : public optimizers::Function {
   
public:

   ScaleFactor();

   ScaleFactor(const optimizers::Function & spectrum,
               double scale_factor=1, bool use_complement=false);

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

   void set_complement_flag(bool use_complement) {
      m_use_complement = use_complement;
      check_complement_usage();
   }

   bool use_complement() const {
      return m_use_complement;
   }

private:
   
   optimizers::Function * m_spectrum;

   bool m_use_complement;

   void init(double scale_factor=1);

   double prefactor() const;

   void check_complement_usage() const;

   void setParRefs();

};

} // namespace Likelihood
