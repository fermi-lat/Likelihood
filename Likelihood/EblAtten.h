/**
 * @brief Decorator class for multiplying underlying spectral models
 * by a user-selectable EBL attenuation.
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/EblAtten.h,v 1.3 2009/06/19 03:39:31 jchiang Exp $
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

   ~EblAtten() throw();

   virtual double value(optimizers::Arg & x) const;
   virtual double derivByParam(optimizers::Arg & x,
                               const std::string & paramName) const;
   virtual Function * clone() const {
      return new EblAtten(*this);
   }

//    /// Set the Parameter value
//    virtual void setParam(const std::string & parName, double value);

   /// Set a Parameter using a Parameter object.
   virtual void setParam(const optimizers::Parameter &param);

//    /// Return the Parameter value by name.
//    virtual double getParamValue(const std::string & parName) const;

//    /// Return the Parameter object by name.
//    virtual const optimizers::Parameter & 
//    getParam(const std::string & parName) const;

//    virtual optimizers::Parameter & parameter(const std::string & parName);

//    virtual optimizers::Parameter & normPar();

//    virtual void getFreeParams(std::vector<optimizers::Parameter> &) const;

//    virtual std::vector<double>::const_iterator
//    setParamValues_(std::vector<double>::const_iterator it);

//    virtual void setParams(const std::vector<optimizers::Parameter> & pars);

//    virtual std::vector<double>::const_iterator
//    setFreeParamValues_(std::vector<double>::const_iterator it);

private:
   
   optimizers::Function * m_spectrum;

   mutable IRB::EblAtten * m_tau;

   void init(double tau_norm=1, double redshift=0, size_t ebl_model=0);

   void setParRefs();

   double attenuation(double energy) const;

};

} // namespace Likelihood
