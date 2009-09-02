/** 
 * @file SkyDirFunction.h
 * @brief Declaration of the SkyDirFunction class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SkyDirFunction.h,v 1.18 2005/02/15 00:34:42 jchiang Exp $
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
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SkyDirFunction.h,v 1.18 2005/02/15 00:34:42 jchiang Exp $
 */
    
class SkyDirFunction : public optimizers::Function {

public:

   SkyDirFunction() {m_init(0., 0.);}
   SkyDirFunction(double ra, double dec) {m_init(ra, dec);}
   SkyDirFunction(const astro::SkyDir &dir);

   const astro::SkyDir & getDir() const {return m_dir;}

   double value(optimizers::Arg &) const {return 0;}

   void setParam(const std::string &paramName, double paramValue) {
      optimizers::Function::setParam(paramName, paramValue);
      update_m_dir(paramName, paramValue);
   }

   double derivByParam(optimizers::Arg &, const std::string &) const
      {return 0;}

   virtual SkyDirFunction * clone() const {
      return new SkyDirFunction(*this);
   }

private:

   void m_init(double ra, double dec);
   
   void update_m_dir(std::string paramName, double paramValue);

   astro::SkyDir::CoordSystem m_coord_type;
   double m_ra;
   double m_dec;

   astro::SkyDir m_dir;

};

} // namespace Likelihood

#endif // Likelihood_SkyDirFunction_h
