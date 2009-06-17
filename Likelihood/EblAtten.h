/**
 * @brief Decorator class for multiplying underlying spectral models
 * by a user-selectable EBL attenuation.
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/EblAtten.h,v 1.1 2009/06/16 05:51:08 jchiang Exp $
 */

namespace Likelihood {

class EblAtten : public optimizers::Function {
   
   EblAtten(const optimizers::Function & spectrum);

   EblAtten(const EblAtten & other);

   ~EblAtten();

   virtual double value(optimizers::Arg & x) const;
   virtual double derivByParam(optimizers::Arg & x,
                               const std::string & paramName) const;
   virtual double integral(optimizers::Arg & xmin, 
                           optimizers::Arg & xmax) const;
   virtual Function * clone() const {
      return new EblAtten(*this);
   }

   /// Set the Parameter value
   virtual void setParam(const std::string & parName, double value);

   /// Set a Parameter using a Parameter object.
   virtual void setParam(const optimizers::Parameter &param);

   /// Return the Parameter value by name.
   virtual double getParamValue(const std::string & parName) const;

   /// Return the Parameter object by name.
   virtual const optimizers::Parameter & 
   getParam(const std::string & parName) const;


private:
   
   optimizers::Function * m_spectrum;

   IRB::EblAtten * m_tau;

   double attenuation(double energy) const;

   bool is_ebl_param(const std::string & paramName) const;

};

} // namespace Likelihood
