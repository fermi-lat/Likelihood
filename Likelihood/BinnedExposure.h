/**
 * @file BinnedExposure.h
 * @brief All-Sky exposure map for use by SourceMap for DiffuseSource 
 * integrations
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef Likelihood_BinnedExposure_h
#define Likelihood_BinnedExposure_h

#include <string>
#include <vector>

#include "astro/SkyDir.h"

#include "map_tools/Exposure.h"

namespace Likelihood {

class BinnedExposure {

public:

   BinnedExposure() {}

   BinnedExposure(const std::vector<double> & energies);

   BinnedExposure(const std::string & filename);

   /// @return Exposure (effective area integrated over time) (cm^2-s)
   /// @param energy True photon energy (MeV)
   /// @param ra Right Ascension of desired sky location (degrees)
   /// @param dec Declination of desired sky location (degrees)
   double operator()(double energy, double ra, double dec) const;

//   void writeOutput(const std::string & filename) const;

private:

   std::vector<double> m_exposureMap;
   std::vector<double> m_ras;
   std::vector<double> m_decs;
   std::vector<double> m_energies;

   unsigned int findIndex(std::vector<double>::const_iterator begin,
                          std::vector<double>::const_iterator end,
                          double value) const;

   void computeMap();

   void linearArray(double xmin, double xmax, unsigned int npts,
                    std::vector<double> &xx) const;

   class Aeff : public map_tools::Exposure::Aeff {
   public:
      Aeff(double energy, int evtType) 
         : m_energy(energy), m_evtType(evtType) {}
      virtual ~Aeff() {}
      virtual double operator()(double cosTheta) const;
   private:
      double m_energy;
      int m_evtType;
      static double s_phi;
   };

};

} // namespace Likelihood

#endif // Likelihood_BinnedExposure_h
