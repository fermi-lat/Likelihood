/**
 * @file ExposureCube.h
 * @brief Exposure time hypercube.
 * @author J. Chiang <jchiang@slacs.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ExposureCube.h,v 1.1 2004/08/23 15:49:51 jchiang Exp $
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ExposureCube.h,v 1.1 2004/08/23 15:49:51 jchiang Exp $
 */

class ExposureCube {

public:

   static ExposureCube * instance() {
      return s_instance;
   }

   static void delete_instance() {
      delete s_instance;
      s_instance = 0;
   }

   static void readExposureCube(std::string filename) {
      s_instance = new ExposureCube();
      facilities::Util::expandEnvVar(&filename);
      s_exposure = new map_tools::Exposure(filename);
   }

   double value(const astro::SkyDir & dir, 
                const map_tools::Exposure::Aeff & aeff) const {
      return (*s_exposure)(dir, aeff);
   }

protected:

   ExposureCube() {}

private:

   static ExposureCube * s_instance;
   static map_tools::Exposure * s_exposure;

};

} // namespace Likelihood

#endif // Likelihood_ExposureCube_h
