/**
 * @file BinnedExposure.h
 * @brief All-Sky exposure map for use by SourceMap for DiffuseSource 
 * integrations
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/BinnedExposure.h,v 1.13 2010/11/24 05:11:26 jchiang Exp $
 */

#ifndef Likelihood_BinnedExposure_h
#define Likelihood_BinnedExposure_h

#include <stdexcept>
#include <string>
#include <vector>

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

class BinnedExposure {

public:

   BinnedExposure();

   BinnedExposure(const CountsMap & cmap, 
                  const Observation & observation, 
                  bool useEbounds=true);

   BinnedExposure(const std::vector<double> & energies,
                  const std::string & proj,
                  const std::string & coordsys,
                  const Observation & observation);

   BinnedExposure(const std::string & filename);

   ~BinnedExposure();

   /// @return Exposure (effective area integrated over time) (cm^2-s)
   /// @param energy True photon energy (MeV)
   /// @param ra Right Ascension of desired sky location (degrees)
   /// @param dec Declination of desired sky location (degrees)
   double operator()(double energy, double ra, double dec) const;

   void writeOutput(const std::string & filename) const;

   const std::vector<double> & energies() const {
      return m_energies;
   }

protected:

// Disable copy constructor and copy assignment operator
   BinnedExposure(const BinnedExposure &) {
      throw std::runtime_error("BinnedExposure copy constructor is disabled");
   }

   BinnedExposure & operator=(const BinnedExposure &) {
      throw std::runtime_error("BinnedExposure copy assignment operator "
                               "is disabled");
      return *this;
   }

   void setMapGeometry(const CountsMap & cmap);

   void setMapGeometry();

private:

   const Observation * m_observation;

   std::vector<float> m_exposureMap;

   std::vector<double> m_energies;

   // Coordinate system parameters to be fed to SkyProj.  These are
   // either set to all-sky values (CAR, CEL) or to match the geometry
   // of the input counts map.
   std::string m_proj_name;
   double m_crpix[2];
   double m_crval[2];
   double m_cdelt[2];
   double m_crota2;
   bool m_isGalactic;

   astro::SkyProj * m_proj;

   std::vector<long> m_naxes;

   void computeMap();
};

} // namespace Likelihood

#endif // Likelihood_BinnedExposure_h
