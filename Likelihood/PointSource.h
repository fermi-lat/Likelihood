/** 
 * @file PointSource.h
 * @brief PointSource class declaration
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/PointSource.h,v 1.45 2004/12/30 00:28:22 jchiang Exp $
 */

#ifndef Likelihood_PointSource_h
#define Likelihood_PointSource_h

#include <utility>

#include "optimizers/Function.h"
#include "optimizers/dArg.h"

#include "Likelihood/Event.h"
#include "Likelihood/ExposureCube.h"
#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/Source.h"

namespace irfInterface {
   class AcceptanceCone;
}

namespace Likelihood {

/** 
 * @class PointSource
 *
 * @brief A gamma-ray point source, e.g., pulsars, AGNs, X-ray binaries.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/PointSource.h,v 1.45 2004/12/30 00:28:22 jchiang Exp $
 */

class PointSource : public Source {

   friend class ExposureMap;

public:

   /// The default constructor does not compute exposure since 
   /// it assumes that the spacecraft data is not available or
   /// that the ROI cuts have not been specified.
   // A later setDir(ra, dec, true) will force the exposure to
   // be computed.
   PointSource() : m_spectrum(0) {setDir(0., 0., false); m_srcType = "Point";}

   /// This constructor does ask for exposure to be computed and 
   /// therefore *requires* the spacecraft data to be available and the 
   /// ROI cuts to be specified beforehand.
   PointSource(double ra, double dec, bool verbose=false) : m_spectrum(0) 
      {setDir(ra, dec, true, verbose);  m_srcType = "Point";}

   PointSource(const PointSource &rhs);

   virtual ~PointSource() {
      delete m_spectrum;
   }

   /// Returns photons/cm^2-s-sr-MeV having been convolved through
   /// the LAT instrument response
   virtual double fluxDensity(const Event &evt) const {
      return fluxDensity(evt.getEnergy(), evt.zAxis(), evt.xAxis(), 
                         evt.getDir(), evt.getType());
   }

   double fluxDensity(double energy, const astro::SkyDir & zAxis,
                      const astro::SkyDir & xAxis,
                      const astro::SkyDir & dir, int eventType=2) const;

   double fluxDensity(double energy, double time,
                      const astro::SkyDir &dir, int eventType=2) const;

   virtual double fluxDensity(double inclination, double phi, double energy, 
                              const astro::SkyDir & appDir, int evtType) const;

   /// Returns the derivative wrt to the named Parameter
   virtual double fluxDensityDeriv(const Event &evt, 
                                   const std::string &paramName) const {
      return fluxDensityDeriv(evt.getEnergy(), evt.zAxis(), evt.xAxis(), 
                              evt.getDir(),evt.getType(), paramName);
   }

   double fluxDensityDeriv(double energy, const astro::SkyDir & zAxis,
                           const astro::SkyDir & xAxis,
                           const astro::SkyDir &dir, int eventType,
                           const std::string &paramName) const;

   double fluxDensityDeriv(double energy, double time,
                           const astro::SkyDir &dir, int eventType,
                           const std::string &paramName) const;

   virtual double fluxDensityDeriv(double inclination, double phi, 
                                   double energy, const astro::SkyDir & appDir,
                                   int evtType, const std::string & paramName)
      const;

   /// Predicted number of photons given RoiCuts and ScData
   virtual double Npred();

   /// Derivative of Npred wrt named Parameter
   virtual double NpredDeriv(const std::string &paramName);

   /// Predicted number of counts within a given energy range
   virtual double Npred(double emin, double emax);

   /// Set source location using J2000 coordinates
   void setDir(double ra, double dec, bool updateExposure=true, 
               bool verbose=true) {
      m_dir = SkyDirFunction(astro::SkyDir(ra, dec));
      m_functions["Position"] = &m_dir;
      if (updateExposure) computeExposure(verbose);
   }

   /// Set source location using Galactic coordinates
   void setGalDir(double l, double b, bool updateExposure=true,
                  bool verbose=true) {
      m_dir = SkyDirFunction(astro::SkyDir(l, b, astro::SkyDir::GALACTIC));
      m_functions["Position"] = &m_dir;
      if (updateExposure) computeExposure(verbose);
   }

   /// Set source location via SkyDir class
   void setDir(const astro::SkyDir &dir, bool updateExposure = true,
               bool verbose=true) {
      m_dir = SkyDirFunction(dir);
      m_functions["Position"] = &m_dir;
      if (updateExposure) computeExposure(verbose);
   }

   const astro::SkyDir & getDir() const {return m_dir.getDir();}

   /// Angular separation between the source direction and dir in radians
   double getSeparation(const astro::SkyDir &dir) 
      {return dir.SkyDir::difference(m_dir.getDir());}

   /// Set the spectral model (@todo Should check that the Parameter
   /// names do not conflict with "longitude" and "latitude" of m_dir)
   void setSpectrum(optimizers::Function *spectrum) {
      m_spectrum = spectrum->clone();
      m_functions["Spectrum"] = m_spectrum;
   }

   virtual Source *clone() const {
      return new PointSource(*this);
   }

   const std::vector<double> & exposure() const {
      return m_exposure;
   }

   virtual double pixelCounts(double emin, double emax,
                              double wtMin, double wtMax) const;

   virtual double pixelCountsDeriv(double emin, double emax, 
                                   double wtMin, double wtMax,
                                   const std::string & paramName) const;

   static bool overlapInterval(const std::pair<double, double> & interval1,
                               std::pair<double, double> & interval2);

private:

   /// Computes the integrated exposure at the PointSource sky location.
   void computeExposure(bool verbose);

   /// Compute the integrated exposure using the provided 
   /// vector of energy values
   void computeExposure(std::vector<double> &energies,
                        std::vector<double> &exposure,
                        bool verbose);

   /// Use a hypercube computed using map_tools.
   void computeExposureWithHyperCube(std::vector<double> &energies,
                                     std::vector<double> &exposure,
                                     bool verbose);

   /// Compute the exposure using integrated exposure times tabulated
   /// as a function of ra, dec, and cos(inclination).
   void computeExposure();

   /// location on the Celestial sphere 
   SkyDirFunction m_dir;

   /// spectral model
   optimizers::Function *m_spectrum;

   /// flag to indicate that static member data has been computed
   static bool s_haveStaticMembers;

   /// vector of energy values for Npred spectrum quadrature
   static std::vector<double> s_energies;

   /// True photon energies for convolving the spectrum with
   /// the energy dispersion.
   static std::vector<double> s_trueEnergies;

   /// integrated exposure at PointSource sky location
   std::vector<double> m_exposure;

   /// method to create a logrithmically spaced grid given RoiCuts
   static void makeEnergyVector(int nee = 100);

   /// Method to compute effective area for the computeExposure time
   /// integrals (when exposure time hypercubes are not available.)
   double sourceEffArea(double energy, double time) const;

#ifndef SWIG
   /**
    * @class Aeff
    *
    * @brief Functor class to calculate PointSource exposure using a
    * map_tools exposure time hypercube.
    *
    */
   class Aeff : public map_tools::Exposure::Aeff {

   public:

      Aeff(double energy, const astro::SkyDir &srcDir);

      virtual ~Aeff() {}

      virtual double operator()(double cos_theta) const;

   private:

      double m_energy;
      astro::SkyDir m_srcDir;

      static std::vector<irfInterface::AcceptanceCone *> s_cones;
      static double s_emin, s_emax;
   };
#endif

};

} //namespace Likelihood

#endif // Likelihood_PointSource_h
