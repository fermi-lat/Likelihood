/**
 * @file RadialProfile.cxx
 * @brief Radial profile for extended sources using ascii file template.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/RadialProfile.cxx,v 1.1 2010/06/19 05:27:31 jchiang Exp $
 */

#include <cmath>

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include "facilities/Util.h"
#include "astro/SkyDir.h"
#include "st_facilities/Util.h"

#include "Likelihood/SkyDirArg.h"
#include "Likelihood/RadialProfile.h"

namespace Likelihood {

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
      m_center = new astro::SkyDir(getParam("ra").getValue(),
                                   getParam("dec").getValue());
   }
   double offset = dir().difference(*m_center)*180./M_PI;
   size_t indx(0);
   if (offset >= m_theta.back()) { // extrapolate beyond last point
      indx = m_theta.size() - 2;
   } else if (offset > m_theta.front()) {
      indx = (std::upper_bound(m_theta.begin(), m_theta.end(), offset)
              - m_theta.begin() - 1);
   }
   return (offset - m_theta.at(indx))/(m_theta.at(indx+1) - m_theta.at(indx))
      *(m_profile.at(indx+1) - m_profile.at(indx)) + m_profile.at(indx);
}

double RadialProfile::derivByParam(optimizers::Arg & x, 
                                   const std::string & parName) const {
   if (parName == "Normalization") {
      return value(x)/getParam(parName).getValue();
   }
   std::ostringstream message;
   message << "RadialProfile: cannot take derivative wrt "
           << "parameter " << parName;
   throw std::runtime_error(message.str());
}

void RadialProfile::setCenter(double ra, double dec) {
   m_center = new astro::SkyDir(ra, dec);
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

void RadialProfile::readTemplateFile(const std::string & template_file) {
   std::vector<std::string> lines;
   bool removeWindowsCRs;
   st_facilities::Util::readLines(template_file, lines, "#", 
                                  removeWindowsCRs=true);
   m_theta.clear();
   m_profile.clear();
   for (size_t i(0); i < lines.size(); i++) {
      std::vector<std::string> tokens;
      facilities::Util::stringTokenize(lines.at(i), " \t", tokens);
      m_theta.push_back(std::atof(tokens.at(0).c_str()));
      m_profile.push_back(std::atof(tokens.at(1).c_str()));
      if (i >= 1 && m_theta.back() < m_theta.at(i-1)) {
         throw std::runtime_error("RadialProfile: angular abscissa points "
                                  "in the template file must be in increasing "
                                  "order.");
      }
   }
}

} // namespace Likelihood
