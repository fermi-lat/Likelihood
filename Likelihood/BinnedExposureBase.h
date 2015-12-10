/**
 * @file BinnedExposureBase.h
 * @brief All-Sky exposure map for use by SourceMap for DiffuseSource 
 * integrations
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/echarles/healpix_changes/Likelihood/Likelihood/BinnedExposureBase.h,v 1.2 2015/03/03 05:59:55 echarles Exp $
 */

#ifndef Likelihood_BinnedExposureBase_h
#define Likelihood_BinnedExposureBase_h

#include <stdexcept>
#include <string>
#include <vector>

#include "Likelihood/ExposureCube.h"

namespace st_app {
   class AppParGroup;
}

namespace astro {
   class SkyProj;
}

namespace Likelihood {

   class CountsMap;
   class Observation;

/**
 * @class BinnedExposureBase
 * @brief This class encapsulates the calculation of and access to 
 * the integral of the effective area over live time.
 *
 * @author J. Chiang
 */

class BinnedExposureBase {

public:

   static double fracDiff(double target, double result) {
     return std::fabs((target - result)/target);
   }
  
   static std::vector<double>::const_iterator
     findNearest(const std::vector<double> & xx, double x, double tol=1e-5);

   BinnedExposureBase();

   BinnedExposureBase(const Observation & observation, 
		      bool useEbounds=true,
		      const st_app::AppParGroup * pars=0);

   BinnedExposureBase(const std::vector<double> & energies,
		      const Observation & observation,
		      const st_app::AppParGroup * pars=0);
   
   BinnedExposureBase(const std::string & filename);

   virtual ~BinnedExposureBase();

   /// @return Exposure (effective area integrated over time) (cm^2-s)
   /// @param energy True photon energy (MeV)
   /// @param ra Right Ascension of desired sky location (degrees)
   /// @param dec Declination of desired sky location (degrees)
   virtual double operator()(double energy, double ra, double dec) const = 0;

   virtual void writeOutput(const std::string & filename) const = 0;

   const std::vector<double> & energies() const {
      return m_energies;
   }

   void setBoundaryFlag(bool enforce_boundaries) {
      m_enforce_boundaries = enforce_boundaries;
   }

   inline bool allSky() const { return m_allSky; }

protected:

// Disable copy constructor and copy assignment operator
   BinnedExposureBase(const BinnedExposureBase &) {
      throw std::runtime_error("BinnedExposureBase copy constructor is disabled");
   }

   BinnedExposureBase & operator=(const BinnedExposureBase &) {
      throw std::runtime_error("BinnedExposureBase copy assignment operator "
                               "is disabled");
      return *this;
   }

//private:

   const Observation * m_observation;

   std::vector<double> m_energies;

   // Coordinate system parameters to be fed to SkyProj.  These are
   // either set to all-sky values (CAR, CEL) or to match the geometry
   // of the input counts map.
   std::string m_proj_name;
   bool m_isGalactic;

   astro::ProjBase * m_proj;

   double m_costhmin;
   double m_costhmax;

   bool m_allSky;
   bool m_enforce_boundaries;

   void setCosThetaBounds(const st_app::AppParGroup & pars);

   class Aeff : public ExposureCube::Aeff {
   public:
      Aeff(double energy, int evtType, const Observation & observation,
           double costhmin, double costhmax) 
         : ExposureCube::Aeff(energy, evtType, observation),
           m_costhmin(costhmin), m_costhmax(costhmax) {}
      virtual double value(double cosTheta, double phi=0) const;
   private:
      double m_costhmin;
      double m_costhmax;
   };
   
};

} // namespace Likelihood

#endif // Likelihood_BinnedExposure_h
