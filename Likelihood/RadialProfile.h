/** 
 * @file RadialProfile.h
 * @brief Radial profile for extended sources using ascii file template.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/RadialProfile.h,v 1.1 2010/06/19 15:17:20 jchiang Exp $
 */

#ifndef Likelihood_RadialProfile_h
#define Likelihood_RadialProfile_h

#include <string>
#include <vector>

namespace astro {
   class SkyDir;
}

#include "optimizers/Function.h"

namespace Likelihood {

/** 
 * @class RadialProfile
 *
 */
    
class RadialProfile : public optimizers::Function {

public:

   RadialProfile();

   RadialProfile(const std::string & template_file);

   RadialProfile(const RadialProfile &);

   virtual ~RadialProfile();

   virtual RadialProfile & operator=(const RadialProfile & rhs);

   double value(optimizers::Arg & dir) const;

   double derivByParam(optimizers::Arg & dir,
                       const std::string & parName) const;

   void readTemplateFile(const std::string & template_file);

   void setCenter(double ra, double dec);

   virtual optimizers::Function * clone() const {
      return new RadialProfile(*this);
   }

private:

   // disable this
   double integral(optimizers::Arg &, optimizers::Arg &) const {return 0;}

   mutable astro::SkyDir * m_center;

   std::vector<double> m_theta;

   std::vector<double> m_profile;

   void init();

};

} // namespace Likelihood

#endif // Likelihood_RadialProfile_h
