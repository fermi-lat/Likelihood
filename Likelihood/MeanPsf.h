/**
 * @file MeanPsf.h
 * @brief Position-dependent Psf averaged over an observation period.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/MeanPsf.h,v 1.15 2010/11/28 23:01:51 jchiang Exp $
 */

#ifndef Likelihood_MeanPsf_h
#define Likelihood_MeanPsf_h

#include "astro/SkyDir.h"

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

   /// @return The value of the psf.
   /// @param energy True photon energy (MeV)
   /// @param theta Angular distance from true source direction (degrees)
   /// @param phi Azimuthal angle about true source direction, 
   ///        currently unused (degrees)
   double operator()(double energy, double theta, double phi=0) const;

   void write(const std::string & filename) const;

   /// Energies (MeV) used for internal representation of psf and
   /// for exposure calculations.
   const std::vector<double> & energies() const {
      return m_energies;
   }

   /// Energy-dependent exposure (cm^2-s) at the selected sky location.
   const std::vector<double> & exposure() const {
      return m_exposure;
   }

   /// @return Exposure at the selected sky location as a function of 
   ///         energy in units of cm^2-s.
   /// @param energy True photon energy (MeV).
   double exposure(double energy) const;

   /// @return Integral over solid angle for an acceptance cone centered
   ///         on the psf.
   /// @param angle Acceptance cone angle (degrees).
   /// @param energy Energy at which to evaluate the integral (MeV).
   double integral(double angle, double energy) const;

   /// @brief Compute the an image of the psf on a longitude-latitude
   ///        grid, assuming the psf center is at (lon0, lat0)
   /// @param energy True energy at which the psf is evaluated (MeV)
   /// @param lon0 Reference longitude (degrees)
   /// @param lat0 Reference latitude (degrees)
   /// @param lons Input longitude array (degrees)
   /// @param lons Input latitude array (degrees)
   /// @param image The output psf image (sr^-1)
   void getImage(double energy, double lon0, double lat0,
                 const std::vector<double> lons,
                 const std::vector<double> lats,
                 std::vector< std::vector<double> > & image) const;

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

   void computeExposure();

   class Psf : public ExposureCube::Aeff {
   public:
      Psf(double separation, double energy, int evtType,
          const Observation & observation) 
         : ExposureCube::Aeff(energy, evtType, observation),
           m_separation(separation) {}
      virtual ~Psf() {}
      virtual double value(double cosTheta, double phi=0) const;
   private:
      double m_separation;
   };

};

}

#endif // Likelihood_MeanPsf_h
