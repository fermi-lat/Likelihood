/**
 * @file BinnedExposure.h
 * @brief All-Sky exposure map for use by SourceMap for DiffuseSource 
 * integrations
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/BinnedExposure.h,v 1.24 2015/12/10 00:57:58 echarles Exp $
 */

#ifndef Likelihood_BinnedExposure_h
#define Likelihood_BinnedExposure_h

#include <stdexcept>
#include <string>
#include <vector>

#include "Likelihood/ExposureCube.h"

// EAC, make a base class for BinnedExposure
#include "Likelihood/BinnedExposureBase.h"

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
 * @class BinnedExposure
 * @brief This class encapsulates the calculation of and access to 
 * the integral of the effective area over live time.
 *
 * @author J. Chiang
 */

class BinnedExposure : public BinnedExposureBase {

public:

   BinnedExposure();

   BinnedExposure(const CountsMap & cmap, 
                  const Observation & observation, 
                  bool useEbounds=true,
                  const st_app::AppParGroup * pars=0);

   BinnedExposure(const std::vector<double> & energies,
                  const Observation & observation,
                  const st_app::AppParGroup * pars=0);

   BinnedExposure(const std::string & filename);

   virtual ~BinnedExposure();

   /// @return Exposure (effective area integrated over time) (cm^2-s)
   /// @param energy True photon energy (MeV)
   /// @param ra Right Ascension of desired sky location (degrees)
   /// @param dec Declination of desired sky location (degrees)
   virtual double operator()(double energy, double ra, double dec) const;

   virtual void get_exposures_for_dir(const astro::SkyDir& dir, 
				      const std::vector<double>& energies, 
				      std::vector<double>& exposures) const;

   virtual void writeOutput(const std::string & filename) const;


protected:

   void setMapGeometry(const CountsMap & cmap, int edisp_bins = 0);

   void setMapGeometry(const st_app::AppParGroup & pars);

   void setMapGeometry();

//private:

   std::vector<float> m_exposureMap;

   // Coordinate system parameters to be fed to SkyProj.  These are
   // either set to all-sky values (CAR, CEL) or to match the geometry
   // of the input counts map.
   double m_crpix[2];
   double m_crval[2];
   double m_cdelt[2];
   double m_crota2;

   std::vector<long> m_naxes;

   void computeMap();

};

} // namespace Likelihood

#endif // Likelihood_BinnedExposure_h
