/** 
 * @file SkyDirFunction.cxx
 * @brief Implementation of SkyDirFunction, a class that encapsulates sky
 * location information in a Function context.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SkyDirFunction.cxx,v 1.11 2003/08/06 20:52:07 jchiang Exp $
 */

#include "Likelihood/SkyDirFunction.h"

namespace Likelihood {

SkyDirFunction::SkyDirFunction(const astro::SkyDir &dir) :
   m_lon(dir.ra()), m_lat(dir.dec()), m_dir(dir) {
   m_maxNumParams = 2;
   m_genericName = "SkyDirFunction";
   addParam("longitude", m_lon, false);
   addParam("latitude", m_lat, false);
}   

void SkyDirFunction::m_init(double lon, double lat) {
   m_lon = lon; 
   m_lat = lat; 
   m_maxNumParams = 2;
   m_genericName = "SkyDirFunction";

   m_dir = astro::SkyDir(lon, lat);

   // add these as fixed parameters 
   // NB: as usual, the specific ordering of Parameters is assumed throughout
   addParam("longitude", m_lon, false);
   addParam("latitude", m_lat, false);
}

void SkyDirFunction::update_m_dir(const std::string paramName, 
                                  double paramValue) 
   throw(optimizers::ParameterNotFound) {
   if (paramName == "longitude") {
      m_lon = paramValue;
   } else if (paramName == "latitude") {
      m_lat = paramValue;
   } else {
      throw optimizers::ParameterNotFound(paramName, getName(), 
                                          "SkyDirFunction::update_m_dir");
   }
   m_dir = astro::SkyDir(m_lon, m_lat);
}

} // namespace Likelihood
