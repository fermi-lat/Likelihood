/**
 * @file BinnedExposure.h
 * @brief All-Sky exposure map for use by SourceMap for DiffuseSource 
 * integrations
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/BinnedExposure.h,v 1.9 2005/11/17 01:47:24 jchiang Exp $
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

private:

   const Observation * m_observation;

   std::vector<float> m_exposureMap;

   std::vector<double> m_energies;

   astro::SkyProj * m_proj;

   std::vector<long> m_naxes;

   void computeMap();

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
