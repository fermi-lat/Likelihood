/**
 * @file ExposureCube.h
 * @brief Exposure time hypercube.
 * @author J. Chiang <jchiang@slacs.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ExposureCube.h,v 1.2 2005/02/27 06:42:24 jchiang Exp $
 */

#ifndef Likelihood_ExposureCube_h
#define Likelihood_ExposureCube_h

#include "facilities/Util.h"

#include "astro/SkyDir.h"

#include "map_tools/Exposure.h"

namespace Likelihood {

/**
 * @class ExposureCube
 * @brief Exposure time as a function of sky position and inclination
 *        wrt the instrument z-axis
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ExposureCube.h,v 1.2 2005/02/27 06:42:24 jchiang Exp $
 */

class ExposureCube {

public:

   ExposureCube() : m_exposure(0), m_haveFile(false) {}

   ~ExposureCube() {
      delete m_exposure;
   }

   void readExposureCube(std::string filename) {
      facilities::Util::expandEnvVar(&filename);
      m_exposure = new map_tools::Exposure(filename);
      m_haveFile = true;
   }

   double value(const astro::SkyDir & dir, 
                const map_tools::Exposure::Aeff & aeff) const {
      return (*m_exposure)(dir, aeff);
   }

   bool haveFile() const {
      return m_haveFile;
   }

private:

   map_tools::Exposure * m_exposure;

   bool m_haveFile;

};

} // namespace Likelihood

#endif // Likelihood_ExposureCube_h
