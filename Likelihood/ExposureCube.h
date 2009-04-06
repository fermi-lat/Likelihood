/**
 * @file ExposureCube.h
 * @brief Exposure time hypercube.
 * @author J. Chiang <jchiang@slacs.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ExposureCube.h,v 1.7 2009/03/16 23:19:06 jchiang Exp $
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ExposureCube.h,v 1.7 2009/03/16 23:19:06 jchiang Exp $
 */

class ExposureCube {

public:

   ExposureCube() : m_exposure(0), m_haveFile(false), m_fileName(""),
                    m_hasPhiDependence(false) {}

   ~ExposureCube() {
      delete m_exposure;
   }

   void readExposureCube(std::string filename) {
      facilities::Util::expandEnvVar(&filename);
      m_fileName = filename;
      m_exposure = new map_tools::Exposure(filename);
      m_haveFile = true;
      m_hasPhiDependence = phiDependence(filename);
   }

   template<class T>
   double value(const astro::SkyDir & dir, const T & aeff) const {
      if (m_hasPhiDependence) {
         AeffWrapper<T> myAeff(aeff);
         return m_exposure->integral(dir, myAeff);
      }
      return (*m_exposure)(dir, aeff);
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

   bool m_hasPhiDependence;

   bool phiDependence(const std::string & filename) const;

   // healpix::CosineBinner::integral assumes that phi is in radians instead of
   // degrees, contrary to established conventions.  This class wraps the 
   // integral(costh, phi) method to do the conversion.
   template <class Aeff>
   class AeffWrapper {
   public:
      AeffWrapper::AeffWrapper(const Aeff & aeff) : m_aeff(aeff) {}
      double integral(double costh, double phi) const {
         phi *= 180./M_PI;
         return m_aeff.integral(costh, phi);
      }
   private:
      const Aeff & m_aeff;
   };

};

} // namespace Likelihood

#endif // Likelihood_ExposureCube_h
