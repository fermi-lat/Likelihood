/**
 * @brief Decorator class for multiplying underlying spectral models
 * by a user-selectable EBL attenuation.
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header$
 */

namespace Likelihood {

class EblAtten : public optimizers::Function {
   
   EblAtten(const optimizers::Function * spectrum);

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

private:
   
   optimizers::Function * m_spectrum;

};

} // namespace Likelihood
