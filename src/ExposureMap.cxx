/**
 * @file ExposureMap.cxx
 * @brief Implementation for the Singleton ExposureMap class. This
 * class encapsulates exposure map information and makes it available
 * for use (primarily) by the DiffuseSource class.
 * @author J. Chiang
 *
 * $Header$
 */

#include "Likelihood/ExposureMap.h"

namespace Likelihood {

ExposureMap * ExposureMap::s_instance = 0;

std::vector<double> ExposureMap::s_ra;
std::vector<double> ExposureMap::s_dec;
std::vector<double> ExposureMap::s_exposure;

void ExposureMap::readExposureFile(std::string exposureFile) {

//     if (!read_FITS_image(exposureFile, ra, dec, ee, expmap)) {
//        computeExposureMap();
//     } else {
//        unpackExposureMap();
//     }
}

ExposureMap * ExposureMap::instance() {
   if (s_instance == 0) {
      s_instance = new ExposureMap();
   }
   return s_instance;
}

} // namespace Likelihood
