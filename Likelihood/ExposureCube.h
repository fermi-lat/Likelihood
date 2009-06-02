/**
 * @file ExposureCube.h
 * @brief Exposure time hypercube.
 * @author J. Chiang <jchiang@slacs.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ExposureCube.h,v 1.14 2009/06/02 16:43:41 jchiang Exp $
 */

#ifndef Likelihood_ExposureCube_h
#define Likelihood_ExposureCube_h

#include "facilities/Util.h"

#include "astro/SkyDir.h"

#include "irfInterface/EfficiencyFactor.h"

#include "map_tools/Exposure.h"

namespace Likelihood {

/**
 * @class ExposureCube
 * @brief Exposure time as a function of sky position and inclination
 *        wrt the instrument z-axis
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ExposureCube.h,v 1.14 2009/06/02 16:43:41 jchiang Exp $
 */

class ExposureCube {

public:

   ExposureCube() : m_exposure(0), m_weightedExposure(0), 
                    m_haveFile(false), m_fileName(""),
                    m_hasPhiDependence(false) {}

   ExposureCube(const ExposureCube & other);

   ~ExposureCube() {
      delete m_exposure;
      delete m_weightedExposure;
   }

   void readExposureCube(std::string filename);

#ifndef SWIG
   template<class T>
   double value(const astro::SkyDir & dir, const T & aeff,
                bool weighted_lt=false) const {
      if (m_hasPhiDependence) {
         AeffWrapper<T> myAeff(aeff);
         if (weighted_lt && m_weightedExposure) {
            return m_weightedExposure->integral(dir, myAeff);
         }
         return m_exposure->integral(dir, myAeff);
      }
      if (weighted_lt && m_weightedExposure) {
         return m_weightedExposure->operator()(dir, aeff);
      }
      return (*m_exposure)(dir, aeff);
   }

   // Compute the exposure with trigger rate- and energy-dependent
   // efficiency corrections.
   template<class T>
   double value(const astro::SkyDir & dir, const T & aeff, 
                double energy) const {
      double factor1, factor2;
      m_efficiencyFactor.getLivetimeFactors(energy, factor1, factor2);
      double exposure(factor1*value(dir, aeff));
      if (factor2 != 0) {
         exposure += factor2*value(dir, aeff, true);
      }
      return exposure;
   }
#endif

   bool haveFile() const {
      return m_haveFile;
   }

   const std::string & fileName() const {
      return m_fileName;
   }

   bool hasPhiDependence() const {
      return m_hasPhiDependence;
   }

private:

#ifndef SWIG
// healpix::CosineBinner::integral assumes that phi is in radians instead of
// degrees, contrary to established conventions.  This class wraps the 
// integral(costh, phi) method to do the conversion.
   template <class Aeff>
   class AeffWrapper {
   public:
      AeffWrapper(const Aeff & aeff) : m_aeff(aeff) {}
      double integral(double costh, double phi) const {
         phi *= 180./M_PI;
         return m_aeff.integral(costh, phi);
      }
   private:
      const Aeff & m_aeff;
   };
#endif

   map_tools::Exposure * m_exposure;
   map_tools::Exposure * m_weightedExposure;

   irfInterface::EfficiencyFactor m_efficiencyFactor;

   bool m_haveFile;

   std::string m_fileName;

   bool m_hasPhiDependence;

   bool phiDependence(const std::string & filename) const;

};

} // namespace Likelihood

#endif // Likelihood_ExposureCube_h
