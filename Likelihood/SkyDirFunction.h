/** 
 * @file SkyDirFunction.h
 * @brief Declaration of the SkyDirFunction class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SkyDirFunction.h,v 1.20 2015/03/03 18:05:36 jchiang Exp $
 */
#ifndef Likelihood_SkyDirFunction_h
#define Likelihood_SkyDirFunction_h

#include "astro/SkyDir.h"

#include "optimizers/Arg.h"
#include "optimizers/Function.h"

namespace Likelihood {

/** 
 * @class SkyDirFunction
 *
 * @brief A class that encapsulates sky location information in a
 * Function context.  This allows Celestial coordinates, i.e., (ra,
 * dec), to be treated by the Likelihood::SourceModel class as
 * Function Parameters.
 *
 */
    
class SkyDirFunction : public optimizers::Function {

public:

   SkyDirFunction(double ra=0, double dec=0);

   SkyDirFunction(const astro::SkyDir & dir);

   const astro::SkyDir & getDir() const {
      return m_dir;
   }

   void setParam(const std::string &paramName, double paramValue) {
      optimizers::Function::setParam(paramName, paramValue);
      update_m_dir(paramName, paramValue);
   }

   virtual Function * clone() const {
      return new SkyDirFunction(*this);
   }

protected:

   virtual double value(const optimizers::Arg &) const {
      return 0;
   }

   virtual double derivByParamImp(const optimizers::Arg &,
                                  const std::string &) const {
      return 0;
   }

private:

   double m_ra;
   double m_dec;

   astro::SkyDir m_dir;

   void init();
   
   void update_m_dir(std::string paramName, double paramValue);
};

} // namespace Likelihood

#endif // Likelihood_SkyDirFunction_h
