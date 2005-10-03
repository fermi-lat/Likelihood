/**
 * @file BinnedExposure.h
 * @brief All-Sky exposure map for use by SourceMap for DiffuseSource 
 * integrations
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/BinnedExposure.h,v 1.6 2005/05/23 19:12:51 jchiang Exp $
 */

#ifndef Likelihood_BinnedExposure_h
#define Likelihood_BinnedExposure_h

#include <string>
#include <vector>

#include "astro/SkyDir.h"

namespace Likelihood {

   class EquinoxRotation;
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

   BinnedExposure(const std::vector<double> & energies,
                  const Observation & observation);

   BinnedExposure(const std::string & filename);

   /// @return Exposure (effective area integrated over time) (cm^2-s)
   /// @param energy True photon energy (MeV)
   /// @param ra Right Ascension of desired sky location (degrees)
   /// @param dec Declination of desired sky location (degrees)
   double operator()(double energy, double ra, double dec) const;

   void writeOutput(const std::string & filename) const;

   void getRotatedImage(double energy, 
                        const std::vector<double> & lons, 
                        const std::vector<double> & lats, 
                        const EquinoxRotation & rot,
                        std::vector< std::vector<double> > & image) const;

private:

   const Observation * m_observation;

   std::vector<float> m_exposureMap;

   std::vector<double> m_ras;
   std::vector<double> m_decs;
   std::vector<double> m_energies;

   unsigned int findIndex(std::vector<double>::const_iterator begin,
                          std::vector<double>::const_iterator end,
                          double value) const;

   void computeMap();

   void linearArray(double xmin, double xmax, unsigned int npts,
                    std::vector<double> &xx) const;

   void fitsReportError(FILE * stream, int status) const;

   class Aeff {
   public:
      Aeff(double energy, int evtType, const Observation & observation) 
         : m_energy(energy), m_evtType(evtType), m_observation(observation) {}
      virtual ~Aeff() {}
      virtual double operator()(double cosTheta) const;
   private:
      double m_energy;
      int m_evtType;
      const Observation & m_observation;
      static double s_phi;
   };

};

} // namespace Likelihood

#endif // Likelihood_BinnedExposure_h
