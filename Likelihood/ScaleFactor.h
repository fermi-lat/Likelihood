/**
 * @brief Decorator class for multiplying underlying spectral models
 * by an overall scale factor (on top of the usual normalization
 * parameter.)  Requested by DM group.
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ScaleFactor.h,v 1.3 2012/04/19 23:39:03 jchiang Exp $
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

   optimizers::Function * spectrum() {
      return m_spectrum;
   }

protected:
   virtual double value(optimizers::Arg & x) const;

   virtual double derivByParamImp(optimizers::Arg & x,
                                  const std::string & paramName) const;

private:
   
   optimizers::Function * m_spectrum;

   bool m_use_complement;

   void init(double scale_factor=1);

   double prefactor() const;

   void check_complement_usage() const;

   void setParRefs();

};

} // namespace Likelihood
