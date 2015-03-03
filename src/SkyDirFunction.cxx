/** 
 * @file SkyDirFunction.cxx
 * @brief Implementation of SkyDirFunction, a class that encapsulates sky
 * location information in a Function context.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SkyDirFunction.cxx,v 1.17 2005/02/02 19:20:03 jchiang Exp $
 */

#include "optimizers/ParameterNotFound.h"

#include "Likelihood/SkyDirFunction.h"

namespace Likelihood {

SkyDirFunction::SkyDirFunction(double ra, double dec)
   : optimizers::Function("SkyDirFunction", 2, ""),
     m_ra(ra), m_dec(dec), m_dir(astro::SkyDir(ra, dec)) {
   init();
}
   
SkyDirFunction::SkyDirFunction(const astro::SkyDir & dir) 
   : optimizers::Function("SkyDirFunction", 2, ""),
     m_ra(dir.ra()), m_dec(dir.dec()), m_dir(dir) {
   init();
}

void SkyDirFunction::init() {
   setName("SkyDirFunction");
   // Add RA, DEC as fixed parameters 
   addParam("RA", m_ra, false);
   parameter("RA").setBounds(-360., 360.);
   setParamAlwaysFixed("RA");
   addParam("DEC", m_dec, false);
   parameter("DEC").setBounds(-90., 90.);
   setParamAlwaysFixed("DEC");
}

void SkyDirFunction::update_m_dir(const std::string paramName, 
                                  double paramValue) {
   if (paramName == "RA") {
      m_ra = paramValue;
   } else if (paramName == "DEC") {
      m_dec = paramValue;
   } else {
      throw optimizers::ParameterNotFound(paramName, getName(), 
                                          "SkyDirFunction::update_m_dir");
   }
   m_dir = astro::SkyDir(m_ra, m_dec);
}

} // namespace Likelihood
