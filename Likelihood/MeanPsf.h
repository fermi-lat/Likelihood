/**
 * @file MeanPsf.h
 * @brief Position-dependent Psf averaged over an observation period.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/MeanPsf.h,v 1.5 2005/02/28 18:38:46 jchiang Exp $
 */

#ifndef Likelihood_MeanPsf_h
#define Likelihood_MeanPsf_h

#include "astro/SkyDir.h"
#include "map_tools/Exposure.h"

#include "Likelihood/Observation.h"

namespace Likelihood {

   class Observation;

/**
 * @class MeanPsf
 *
 */

class MeanPsf {

public:

   MeanPsf(double ra, double dec, const std::vector<double> & energies,
           const Observation & observation) 
      : m_srcDir(ra, dec, astro::SkyDir::EQUATORIAL),
        m_energies(energies), m_observation(observation) {
      init();
   }

   MeanPsf(const astro::SkyDir & srcDir, const std::vector<double> & energies,
           const Observation & observation) 
      : m_srcDir(srcDir), m_energies(energies), m_observation(observation) {
      init();
   }

   MeanPsf() : m_observation(Observation()) {}

   /// @return The value of the psf.
   /// @param energy True photon energy (MeV)
   /// @param theta Angular distance from true source direction (degrees)
   /// @param phi Azimuthal angle about true source direction, 
   ///        currently unused (degrees)
   double operator()(double energy, double theta, double phi=0) const;

   void write(const std::string & filename) const;

   /// Energy-dependent exposure (cm^2-s) at the selected sky location.
   const std::vector<double> & exposure() const {
      return m_exposure;
   }

private:

   static std::vector<double> s_separations;

   astro::SkyDir m_srcDir;

   std::vector<double> m_energies;

   const Observation & m_observation;

   std::vector<double> m_psfValues;

   std::vector<double> m_exposure;

   void init();

   void createLogArray(double xmin, double xmax, unsigned int npts,
                       std::vector<double> & xx) const;


   class Psf : public map_tools::Exposure::Aeff {
   public:
      Psf(double separation, double energy, int evtType,
          const Observation & observation) 
         : m_separation(separation), m_energy(energy), m_evtType(evtType),
           m_observation(observation) {}
      virtual ~Psf() {}
      virtual double operator()(double cosTheta) const;
   private:
      double m_separation;
      double m_energy;
      int m_evtType;
      const Observation & m_observation;
      static double s_phi;
   };

   class Aeff : public map_tools::Exposure::Aeff {
   public:
      Aeff(double energy, int evtType, const Observation & observation) 
         : m_energy(energy), m_evtType(evtType), m_observation(observation) {}
      virtual ~Aeff() {}
      virtual double operator()(double cosTheta) const;
   private:
      double m_separation;
      double m_energy;
      int m_evtType;
      const Observation & m_observation;
      static double s_phi;
   };

};

}

#endif // Likelihood_MeanPsf_h
