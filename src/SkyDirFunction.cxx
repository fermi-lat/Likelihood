/** 
 * @file SkyDirFunction.cxx
 * @brief Implementation of SkyDirFunction, a class that encapsulates sky
 * location information in a Function context.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SkyDirFunction.cxx,v 1.15 2004/12/22 06:06:48 jchiang Exp $
 */

#include "optimizers/ParameterNotFound.h"

#include "Likelihood/SkyDirFunction.h"

namespace Likelihood {

SkyDirFunction::SkyDirFunction(const astro::SkyDir &dir) :
   m_ra(dir.ra()), m_dec(dir.dec()), m_dir(dir) {
   m_maxNumParams = 2;
   m_genericName = "SkyDirFunction";
   m_functionName = "SkyDirFunction";
   addParam("RA", m_ra, false);
   parameter("RA").setBounds(-360., 360.);
   addParam("DEC", m_dec, false);
   parameter("DEC").setBounds(-90., 90.);
}   

void SkyDirFunction::m_init(double ra, double dec) {
   m_ra = ra; 
   m_dec = dec; 
   m_maxNumParams = 2;
   m_genericName = "SkyDirFunction";
   m_functionName = "SkyDirFunction";

   m_dir = astro::SkyDir(ra, dec);

   // add these as fixed parameters 
   // NB: as usual, the specific ordering of Parameters is assumed throughout
   addParam("RA", m_ra, false);
   parameter("RA").setBounds(-360., 360.);
   addParam("DEC", m_dec, false);
   parameter("DEC").setBounds(-90., 90.);
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
