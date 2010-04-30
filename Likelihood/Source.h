/** 
 * @file Source.h
 * @brief Source base class declaration
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Source.h,v 1.45 2010/03/07 18:23:57 jchiang Exp $
 */

#ifndef Likelihood_Source_h
#define Likelihood_Source_h

#include <iostream>
#include <map>
#include <stdexcept>

#include "optimizers/dArg.h"
#include "optimizers/Function.h"

#include "Likelihood/Event.h"

namespace astro {
   class SkyDir;
}

namespace Likelihood {

   class Observation;

/** 
 * @class Source
 *
 * @brief Base class for gamma-ray sources.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Source.h,v 1.45 2010/03/07 18:23:57 jchiang Exp $
 */

class Source {

public:
   typedef std::pair<bool, double> CachedResponse;
    
   Source(const Observation * observation=0);

   Source(const Source &rhs);

   virtual ~Source() {
      delete m_spectrum;
   }

   /// @return photons/cm^2-s-sr-MeV having been convolved through
   /// the LAT instrument response
   virtual double fluxDensity(const Event &evt, 
                              CachedResponse* cResp = 0) const = 0;

   /// @return fluxDensity in instrument coordinates (photons/cm^2-s-sr-MeV)
   /// @param inclination angle of source direction wrt the instrument
   ///        z-axis (degrees)
   /// @param phi azimuthal angle of source direction wrt instrument
   ///        z- and x-axes (degrees)
   /// @param energy True energy of photon (MeV)
   /// @param appDir Apparent photon direction
   /// @param evtType Event type, i.e., front- vs back-converting event, 
   ///        0 vs 1
   virtual double fluxDensity(double inclination, double phi, double energy, 
                              const astro::SkyDir & appDir, 
                              int evtType, CachedResponse* cResp = 0) const = 0;

   /// Derivatives of fluxDensity wrt model Parameters
   virtual double fluxDensityDeriv(const Event &evt, 
                                   const std::string &paramName,
				   CachedResponse* cResp = 0) const = 0;

   virtual double fluxDensityDeriv(double inclination, double phi, 
                                   double energy, const astro::SkyDir & appDir,
                                   int evtType, 
                                   const std::string & paramName,
				   CachedResponse* cResp = 0) const = 0;

   /// Predicted number of photons.
   virtual double Npred();

   /// Derivative of Npred wrt named Parameter
   virtual double NpredDeriv(const std::string & paramName);

   /// Predicted number of counts within a specified energy range.
   virtual double Npred(double emin, double emax) const;

   /// Set the spectral model (should also check that the Parameter
   /// names do not conflict with "longitude" and "latitude" of m_dir)
   virtual void setSpectrum(optimizers::Function * spectrum) {
      m_spectrum = spectrum->clone();
      m_functions["Spectrum"] = m_spectrum;
   }

   /// Set the spectral model by name.
   virtual void setSpectrum(const std::string & functionName);
                       
   virtual void setName(const std::string & name) {
      m_name = name;
   }

   /// @return true if all spectral parameters are fixed.
   bool fixedSpectrum() const;

   /// Access a unique source identifier.
   virtual const std::string & getName() const {
      return m_name;
   }

   typedef std::map<std::string, optimizers::Function *> FuncMap;

   /// @return A mutable reference to the m_functions map
   virtual FuncMap & getSrcFuncs() {
      return m_functions;
   }

   virtual const FuncMap & getSrcFuncs() const {
      return m_functions;
   }

   virtual Source * clone() const {
      return 0;
   }

   /// The Source type (e.g., Diffuse vs Point)
   virtual const std::string & getType() const {
      return m_srcType;
   }

   /// Integrate the product of the source spectrum with the given
   /// SourceMap pixel values.
   virtual double pixelCounts(double emin, double emax, 
                              double wtMin, double wtMax) const;

   virtual double pixelCountsDeriv(double emin, double emax, 
                                   double wtMin, double wtMax,
                                   const std::string & paramName) const;

   virtual const std::vector<double> & exposure() const = 0;

   virtual const optimizers::Function & spectrum() const = 0;

   virtual double flux() const = 0;

   virtual double fluxDeriv(const std::string & parName) const = 0;

   /// @return Photon flux integrated over the given energy range.
   /// Units are #/cm^2/s
   virtual double flux(double emin, double emax, size_t npts=100) const = 0;

   /// @return Derivative of integrated photon flux wrt the named parameter
   /// over the given energy range.
   virtual double fluxDeriv(const std::string & parName, double emin,
                            double emax, size_t npts=100) const = 0;

   /// @return Energy flux integrated over the ROI energy bounds. 
   /// Units are MeV/cm^2/s
   virtual double energyFlux() const = 0;

   /// @return Derivative of integrated energy flux wrt the named parameter
   virtual double energyFluxDeriv(const std::string & parName) const = 0;

   /// @return Energy flux integrated over the given energy range.
   /// Units are MeV/cm^2/s
   virtual double energyFlux(double emin, double emax,
                             size_t npts=100) const = 0;

   /// @return Derivative of integrated energy flux wrt the named parameter
   /// over the given energy range.
   virtual double energyFluxDeriv(const std::string & parName, double emin,
                                  double emax, size_t npts=100) const = 0;

   const Observation * observation() const {
      return m_observation;
   }

   void setObservation(const Observation * observation) {
      m_observation = observation;
   }

protected:

   /// A unique source name.
   std::string m_name;

   /// Source type --- "Diffuse" or "Point".
   std::string m_srcType;

   /// Flag to indicate if energy dispersion is to be used.
   bool m_useEdisp;

   /// Map of Functions describing this source.
   std::map<std::string, optimizers::Function *> m_functions;

   /// spectral model
   optimizers::Function * m_spectrum;

   const Observation * m_observation;

   /// Angle integrated diffuse exposure as a function of
   /// RoiCuts::energies()
   std::vector<double> m_exposure;

   static double powerlaw_integral_est(double x1, double x2, 
                                       double y1, double y2, 
                                       double wt1, double wt2);

/// Nested classes for computing photon flux derivatives and energy
/// fluxes and flux derivatives.

/**
 * @class FluxDeriv
 * @brief Functor class that wraps a Function to provide an interface
 * to that function's partial derivative wrt a named parameter.
 */
   class FluxDeriv : public optimizers::Function {
   public:
      FluxDeriv(const optimizers::Function & func, const std::string & parName) 
         : m_func(func), m_parName(parName) {}
      virtual double value(optimizers::Arg & x) const {
         return m_func.derivByParam(x, m_parName);
      }
      virtual double derivByParam(optimizers::Arg &,
                                  const std::string &) const {
         throw std::runtime_error("FluxDeriv::deriveByParam not implemented");
      }
   protected:
      virtual Function * clone() const {
         return 0;
      }
   private:
      const optimizers::Function & m_func;
      std::string m_parName;
   };

/**
 * @class EnergyFlux
 * @brief Functor class to be used to compute energy fluxes.
 */
   class EnergyFlux : public optimizers::Function {
   public:
      EnergyFlux(const optimizers::Function & func) : m_func(func) {}
      virtual double value(optimizers::Arg & x) const {
         double energy(dynamic_cast<optimizers::dArg &>(x).getValue());
         return energy*m_func(x);
      }
      virtual double derivByParam(optimizers::Arg & x,
                                  const std::string & parname) const {
         double energy(dynamic_cast<optimizers::dArg &>(x).getValue());
         return energy*m_func.derivByParam(x, parname);
      }
   protected:
      virtual Function * clone() const {
         return 0;
      }
   private:
      const optimizers::Function & m_func;
   };

/**
 * @class EnergyFluxDeriv
 * @brief Functor class to be used to compute energy flux derivatives
 * wrt fit parameters.
 */
   class EnergyFluxDeriv : public optimizers::Function {
   public:
      EnergyFluxDeriv(const optimizers::Function & func,
                      const std::string & parName) 
         : m_func(func), m_parName(parName) {}
      virtual double value(optimizers::Arg & x) const {
         double energy(dynamic_cast<optimizers::dArg &>(x).getValue());
         return energy*m_func.derivByParam(x, m_parName);
      }
      virtual double derivByParam(optimizers::Arg &,
                                  const std::string &) const {
         throw std::runtime_error("EnergyFluxDeriv::derivByParam: "
                                  "not implemented");
      }
   protected:
      virtual Function * clone() const {
         return 0;
      }
   private:
      const optimizers::Function & m_func;
      std::string m_parName;
   };

};

} // namespace Likelihood

#endif // Likelihood_Source_h
