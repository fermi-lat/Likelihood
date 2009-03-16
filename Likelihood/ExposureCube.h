/**
 * @file ExposureCube.h
 * @brief Exposure time hypercube.
 * @author J. Chiang <jchiang@slacs.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ExposureCube.h,v 1.5 2006/04/24 22:31:07 jchiang Exp $
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ExposureCube.h,v 1.5 2006/04/24 22:31:07 jchiang Exp $
 */

class ExposureCube {

public:

   ExposureCube() : m_exposure(0), m_haveFile(false), m_fileName("") {}

   ~ExposureCube() {
      delete m_exposure;
   }

   void readExposureCube(std::string filename) {
      facilities::Util::expandEnvVar(&filename);
      m_fileName = filename;
      m_exposure = new map_tools::Exposure(filename);
      m_haveFile = true;
   }

   template<class T>
   double value(const astro::SkyDir & dir, const T & aeff) const {
//      return (*m_exposure)(dir, aeff);
      return m_exposure->integral(dir, aeff);
   }

   bool haveFile() const {
      return m_haveFile;
   }

   const std::string & fileName() const {
      return m_fileName;
   }

private:

   map_tools::Exposure * m_exposure;

   bool m_haveFile;

   std::string m_fileName;

};

} // namespace Likelihood

#endif // Likelihood_ExposureCube_h
