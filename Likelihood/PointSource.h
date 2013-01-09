/** 
 * @file PointSource.h
 * @brief PointSource class declaration
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/PointSource.h,v 1.72 2012/06/27 20:31:50 jchiang Exp $
 */

#ifndef Likelihood_PointSource_h
#define Likelihood_PointSource_h

#include <utility>

#include "optimizers/Function.h"
#include "optimizers/dArg.h"

#include "Likelihood/Event.h"
#include "Likelihood/ExposureCube.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/Source.h"

namespace st_stream {
   class StreamFormatter;
}

namespace irfInterface {
   class AcceptanceCone;
}

namespace Likelihood {

   class Observation;
   class ResponseFunctions;
   class RoiCuts;
   class ScData;

/** 
 * @class PointSource
 *
 * @brief A gamma-ray point source, e.g., pulsars, AGNs, X-ray binaries.
 *
 */

class PointSource : public Source {

public:

   /// The default constructor does not compute exposure since 
   /// it assumes that the spacecraft data is not available or
   /// that the ROI cuts have not been specified.
   // A later setDir(ra, dec, true) will force the exposure to
   // be computed.
   PointSource(const Observation * observation=0);

   /// This constructor does ask for exposure to be computed and 
   /// therefore *requires* the spacecraft data to be available and the 
   /// ROI cuts to be specified beforehand.
   PointSource(double ra, double dec, const Observation & observation,
               bool verbose=false);

   PointSource(const PointSource & rhs);

   virtual ~PointSource();

   /// Returns photons/cm^2-s-sr-MeV having been convolved through
   /// the LAT instrument response
   virtual double fluxDensity(const Event & evt,
                              CachedResponse* cResp = 0) const {
      return fluxDensity(evt.getEnergy(), evt.zAxis(), evt.xAxis(), 
                         evt.getDir(), evt.getType(), evt.getArrTime(),
                         cResp);
   }

   virtual double fluxDensity(double inclination, double phi, double energy, 
                              const astro::SkyDir & appDir, int evtType,
                              double time, CachedResponse* cResp = 0) const;

   /// Returns the derivative wrt to the named Parameter
   virtual double fluxDensityDeriv(const Event & evt, 
                                   const std::string & paramName,
				   CachedResponse* cResp = 0) const {
      return fluxDensityDeriv(evt.getEnergy(), evt.zAxis(), evt.xAxis(), 
                              evt.getDir(),evt.getType(), evt.getArrTime(),
                              paramName, cResp);
   }

   virtual double fluxDensityDeriv(double inclination, double phi, 
                                   double energy, const astro::SkyDir & appDir,
                                   int evtType, double time,
                                   const std::string & paramName,
				   CachedResponse* cResp = 0) const;

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

   const astro::SkyDir & getDir() const {
      return m_dir.getDir();
   }

   /// Angular separation between the source direction and dir in radians
   double getSeparation(const astro::SkyDir &dir) const {
      return dir.SkyDir::difference(m_dir.getDir());
   }

   virtual Source *clone() const {
      return new PointSource(*this);
   }

   virtual void computeExposure(const std::vector<double> & energies,
                                bool verbose=false);

   /// Compute the integrated exposure using the provided 
   /// vector of energy values
   static void computeExposure(const astro::SkyDir & dir,
                               const std::vector<double> & energies,
                               const Observation & observation, 
                               std::vector<double> & exposure,
                               bool verbose);

   /// Use a hypercube computed using map_tools.
   static void 
   computeExposureWithHyperCube(const astro::SkyDir & dir, 
                                const std::vector<double> & energies,
                                const Observation & observation, 
                                std::vector<double> & exposure,
                                bool verbose);

   virtual const std::vector<double> & exposure() const {
      return m_exposure;
   }

   // virtual const optimizers::Function & spectrum() const {
   //    return * m_spectrum;
   // }

   // virtual optimizers::Function & spectrum() {
   //    return *m_spectrum;
   // }

   /// @return Photon flux integrated over the ROI energy bounds. 
   /// Units are #/cm^2/s
   double flux() const;

   /// @return Derivative of integrated photon flux wrt the named parameter
   double fluxDeriv(const std::string & parName) const;

   /// @return Photon flux integrated over the given energy range.
   /// Units are #/cm^2/s
   double flux(double emin, double emax, size_t npts=100) const;

   /// @return Derivative of integrated photon flux wrt the named parameter
   /// over the given energy range.
   double fluxDeriv(const std::string & parName, 
                    double emin, double emax, size_t npts=100) const;

   /// @return Energy flux integrated over the ROI energy bounds. 
   /// Units are MeV/cm^2/s
   double energyFlux() const;

   /// @return Derivative of integrated energy flux wrt the named parameter
   double energyFluxDeriv(const std::string & parName) const;

   /// @return Energy flux integrated over the given energy range.
   /// Units are MeV/cm^2/s
   double energyFlux(double emin, double emax, size_t npts=100) const;

   /// @return Derivative of integrated energy flux wrt the named parameter
   /// over the given energy range.
   double energyFluxDeriv(const std::string & parName, 
                          double emin, double emax, size_t npts=100) const;

private:

   /// location on the Celestial sphere 
   SkyDirFunction m_dir;

   /// True photon energies for convolving the spectrum with
   /// the energy dispersion.
   static std::vector<double> s_trueEnergies;

   /// Computes the integrated exposure at the PointSource sky location.
   void computeExposure(bool verbose);

   double fluxDensity(double energy, const astro::SkyDir & zAxis,
                      const astro::SkyDir & xAxis,
                      const astro::SkyDir & dir, 
                      int eventType=2,
                      double time=0,
		      CachedResponse* cResp = 0) const;

   double fluxDensity(double energy, double time,
                      const astro::SkyDir &dir, int eventType=2,
		      CachedResponse* cResp = 0) const;

   double fluxDensityDeriv(double energy, const astro::SkyDir & zAxis,
                           const astro::SkyDir & xAxis,
                           const astro::SkyDir &dir, int eventType,
                           double time,
                           const std::string &paramName,
			   CachedResponse* cResp = 0) const;

   double fluxDensityDeriv(double energy, double time,
                           const astro::SkyDir &dir, int eventType,
                           const std::string &paramName,
			   CachedResponse* cResp = 0) const;

   /// method to create a logrithmically spaced grid given RoiCuts
   static void makeEnergyVector(int nee = 100);

   /// Method to compute effective area for the computeExposure time
   /// integrals (when exposure time hypercubes are not available.)
   static double sourceEffArea(const astro::SkyDir & srcDir,
                               double energy, double time,
                               const ScData & scData,
                               const RoiCuts & roiCuts,
                               const ResponseFunctions & respFuncs);

#ifndef SWIG
   /**
    * @class Aeff
    *
    * @brief Functor class to calculate PointSource exposure using a
    * map_tools exposure time hypercube.
    *
    */
   class Aeff : public ExposureCube::AeffBase {

   public:

      Aeff(double energy, const astro::SkyDir &srcDir,
           const RoiCuts & roiCuts, const ResponseFunctions & respFuncs,
           double time=0, bool usePhiDependence=true);
      virtual ~Aeff() {}

   private:

      double m_energy;
      astro::SkyDir m_srcDir;

      const ResponseFunctions & m_respFuncs;

      double m_time;

      bool m_usePhiDependence;

      std::vector<irfInterface::AcceptanceCone *> m_cones;
      double m_emin, m_emax;

      virtual double value(double cos_theta, double phi=0) const;
   };
#endif

};

} //namespace Likelihood

#endif // Likelihood_PointSource_h
