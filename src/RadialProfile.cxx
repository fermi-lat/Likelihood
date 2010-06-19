/**
 * @file RadialProfile.cxx
 * @brief Radial profile for extended sources using ascii file template.
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>

#include "astro/SkyDir.h"
#include "RadialProfile.h"

//namespace Likelihood {

RadialProfile::RadialProfile() : m_center(0) {
   init();
}

RadialProfile::RadialProfile(const std::string & template_file) 
  : m_center(0) {
   init();
   readTemplateFile(template_file);
}

RadialProfile::RadialProfile(const RadialProfile & other) 
   : optimizers::Function(other),
     m_center(new astro::SkyDir(*other.m_center)) ,
     m_theta(other.m_theta), 
     m_profile(other.m_profile) {}

RadialProfile::~RadialProfile() {
   delete m_center;
}

RadialProfile & RadialProfile::operator=(const RadialProfile & rhs) {
   if (this != &rhs) {
      optimizers::Function::operator=(rhs);
      m_center = new astro::SkyDir(*rhs.m_center);
      m_theta = rhs.m_theta;
      m_profile = rhs.m_profile;
   }
   return *this;
}

double RadialProfile::value(optimizers::Arg & x) const {
   SkyDirArg & dir = dynamic_cast<SkyDirArg &>(x);
   if (m_center == 0) {
      m_center = new astro::SkyDir(parameter("ra").getValue(),
                                   parameter("dec").getValue());
   }
   double offset = dir.difference(*m_center)*180./M_PI;
   size_t indx = std::upper_bound(m_theta.begin(), m_theta.end(), offset);
}

void RadialProfile::init() {
   setMaxNumParams(3);

   addParam("Normalization", 1, false);
   parameter("Normalization").setBounds(0, 10);
   setParamAlwaysFixed("Normalization");

   addParam("ra", 0, false);
   parameter("ra").setBounds(-360, 360);
   setParamAlwaysFixed("ra");

   addParam("dec", 0, false);
   parameter("dec").setBounds(-90, 90);
   setParamAlwaysFixed("dec");

   m_genericName = "RadialProfile";
   m_normParName = "Normalization";
   m_funcType = Addend;
   m_argType = "";
}

} // namespace Likelihood
