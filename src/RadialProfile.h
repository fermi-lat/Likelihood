/** 
 * @file RadialProfile.h
 * @brief Radial profile for extended sources using ascii file template.
 * @author J. Chiang
 *
 * $Header$
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

   virtual optimizers::Function * clone() const {
      return new RadialProfile(*this);
   }

private:

   // disable this
   double integral(Arg &, Arg &) const {return 0;}

   astro::SkyDir * m_center;

   std::vector<double> m_theta;

   std::vector<double> m_profile;

   void init();

   void readTemplateFile(const std::string & template_file);

};

} // namespace Likelihood

#endif // Likelihood_RadialProfile_h
