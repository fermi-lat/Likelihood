/**
 * @brief Decorator class for multiplying underlying spectral models
 * by a user-selectable EBL attenuation.
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/EblAtten.h,v 1.6 2015/03/03 18:05:36 jchiang Exp $
 */

#include "optimizers/Function.h"

#include "eblAtten/EblAtten.h"

namespace Likelihood {

class EblAtten : public optimizers::Function {
   
public:

   EblAtten();

   EblAtten(const optimizers::Function & spectrum,
            double tau_norm=1, double redshift=0, 
            size_t ebl_model=0);

   EblAtten(const EblAtten & other);

   EblAtten & operator=(const EblAtten & rhs);

   ~EblAtten() throw();

   virtual Function * clone() const {
      return new EblAtten(*this);
   }

   /// Set a Parameter using a Parameter object.  This version
   /// preserves the references to the m_spectrum parameters.
   virtual void setParam(const optimizers::Parameter &param);

protected:

   virtual double value(const optimizers::Arg & x) const;

   virtual double derivByParamImp(const optimizers::Arg & x,
                                  const std::string & paramName) const;

private:
   
   optimizers::Function * m_spectrum;

   mutable IRB::EblAtten * m_tau;

   void init(double tau_norm=1, double redshift=0, size_t ebl_model=0);

   void setParRefs();

   double attenuation(double energy) const;

};

} // namespace Likelihood
