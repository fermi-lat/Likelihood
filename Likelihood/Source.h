/** 
 * @file Source.h
 * @brief Source base class declaration
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Source.h,v 1.54 2016/09/09 21:11:52 echarles Exp $
 */

#ifndef Likelihood_Source_h
#define Likelihood_Source_h

#include <iostream>
#include <map>
#include <stdexcept>

#include "optimizers/dArg.h"
#include "optimizers/Function.h"

#include "Likelihood/Event.h"
#include "Likelihood/TrapQuad.h"

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
    */

   class Source {
     
   public:
     
     typedef enum { 
       //! Point Source
       Point = 0,
       //! Diffuse sources
       Diffuse = 1,
       //! Composite sources
       Composite = 2,
       //! Unknown
       Unknown = 3 } SourceType;      
     
     /* This is for backwards compatibilty, to convert the enum to a string */
     static const std::string& sourceTypeName(SourceType t);
   
     /* FIXME, what exactly is this */
     typedef std::pair<bool, double> CachedResponse;

     /* Type for the map of the functions need to define a source. */
     typedef std::map<std::string, optimizers::Function *> FuncMap;
    
   public:
     
     /* Default c'tor */
     Source(const Observation * observation=0);

     /* Copy c'tor */
     Source(const Source &rhs);

     /* D'tor */
     virtual ~Source() {
       delete m_spectrum;
     }

     /* Clone function, implemented by sub-classes */
     virtual Source * clone() const {
       return 0;
     }


     /* ------------- Simple access functions ---------------------- */

     /// Access a unique source identifier.
     inline const std::string & getName() const { return m_name; }

     /// @return A mutable reference to the m_functions map
     inline FuncMap & getSrcFuncs() { return m_functions; }
     /// @return A const reference to the m_functions map
     inline const FuncMap & getSrcFuncs() const { return m_functions; }

     /// The Source type (e.g., Diffuse vs Point)
     inline SourceType srcType() const {  return m_srcType; }

     /// This is for backwards compatibility, return a string for the source type
     inline const std::string & getType() const { return sourceTypeName(m_srcType); }

     /// @return A mutable reference to the Spectrum function      
     inline optimizers::Function & spectrum() { return *m_spectrum; }
     /// @return A const reference to the Spectrum function 
     inline const optimizers::Function & spectrum() const { return *m_spectrum; }

     /// @return A pointer to the observation container 
     inline const Observation * observation() const { return m_observation; }

     /// @return Flag to enable energy dispersion
     bool use_edisp() const { return m_useEdisp; }

     /// @return Cached vector of exposures as a function of energy
     inline const std::vector<double> & exposure() const { return m_exposure; }
 

     /* ----------------- Simple setter functions ------------------------- */
   
     /// Set the spectral model
     /// FIXME (should also check that the Parameter
     /// names do not conflict with "longitude" and "latitude" of m_dir)
     virtual void setSpectrum(optimizers::Function * spectrum) {
       m_spectrum = spectrum->clone();
       m_functions["Spectrum"] = m_spectrum;
     }
     
     /// Set the spectral model by name.
     virtual void setSpectrum(const std::string & functionName);
     
     /// Set the name 
     inline void setName(const std::string & name) { m_name = name; }   

     /// Set the observation container
     inline void setObservation(const Observation * observation) { m_observation = observation; }

     /// Set the energy dispersion flag
     inline void set_edisp_flag(bool useEdisp) { m_useEdisp = useEdisp; }


     /* ------------------ Functions for unbinned likelihood --------------- */


     /// @return photons/cm^2-s-sr-MeV having been convolved through
     /// the LAT instrument response for a particular event
     /// @param evt container for event data
     /// @param cResp Cached instrument response
     virtual double fluxDensity(const Event & evt, 
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
     /// @param time Event arrival time
     /// @param cResp Cached instrument response
     virtual double fluxDensity(double inclination, double phi, double energy, 
				const astro::SkyDir & appDir, 
				int evtType, double time, 
				CachedResponse* cResp = 0) const = 0;

     /// @return derivative of fluxDensity wrt a model Parameter
     /// @param evt container for event data
     /// @param paramName name of the parameter in question
     /// @param cResp Cached instrument response   
     virtual double fluxDensityDeriv(const Event & evt, 
				     const std::string & paramName,
				     CachedResponse* cResp = 0) const = 0;
     
     
     /// @return derivative of fluxDensity wrt a model Parameter
     /// @param inclination angle of source direction wrt the instrument
     ///        z-axis (degrees)
     /// @param phi azimuthal angle of source direction wrt instrument
     ///        z- and x-axes (degrees)
     /// @param energy True energy of photon (MeV)
     /// @param appDir Apparent photon direction
     /// @param evtType Event type, i.e., front- vs back-converting event, 
     ///        0 vs 1
     /// @param time Event arrival time
     /// @param paramName name of the parameter in question
     /// @param cResp Cached instrument response
     virtual double fluxDensityDeriv(double inclination, double phi, 
				     double energy, const astro::SkyDir & appDir,
				     int evtType, double time, 
				     const std::string & paramName,
				     CachedResponse* cResp = 0) const = 0;

     

     /// FIXME, what exactly does this do?
     virtual void computeExposure(const std::vector<double> & energies,
				  bool verbose=false) = 0;


     /* ------------------ Functions for binned likelihood --------------- */
     
     /// Integrate the product of the source spectrum with the given
     /// SourceMap pixel values.
     double pixelCounts(double emin, double emax, 
			double wtMin, double wtMax) const;
     
     /// Integrate the derivative of the source spectrum with the given
     /// SourceMap pixel values.
     double pixelCountsDeriv(double emin, double emax, 
			     double wtMin, double wtMax,
			     const std::string & paramName) const;
     
  
     /* ------------ Functions for both binned and unbinned likelihood ------------ */
     
     /// @return true if all spectral parameters are fixed.
     bool fixedSpectrum() const;
      
     /// Predicted number of photons for this source (for the event selection)
     virtual double Npred();

     /// Derivative of Npred wrt named Parameter
     virtual double NpredDeriv(const std::string & paramName);

     /// Predicted number of counts within a specified energy range (for the event selection)
     virtual double Npred(double emin, double emax) const;
     
     
     /// Derivative of Npred within a specified energy range (for the event selection) wrt named Paramter
     virtual double NpredDeriv(const std::string & paramName,
                             double emin, double emax) const;

     /// @return Photon flux integrated over the ROI energy bounds. 
     /// Units are #/cm^2/s
     virtual double flux() const = 0;

     /// @return Derivative of integrated photon flux wrt the named parameter
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
     
     
     /// @return the spatial integral of the source.  
     /// This is here to allow for different normalizations in diffuse maps.
     virtual double angularIntegral(double energy) const { return 1.0; }

   protected:
     
     /// A unique source name.
     std::string m_name;
     
     /// Source type --- "Diffuse" or "Point".
     SourceType m_srcType;
     
     /// Flag to indicate if energy dispersion is to be used.
     bool m_useEdisp;
     
     /// Map of Functions describing this source.
     std::map<std::string, optimizers::Function *> m_functions;
     
     /// spectral model
     optimizers::Function * m_spectrum;
     
     /// Container for infomation about the observation
     const Observation * m_observation;
     
     /// Angle integrated diffuse exposure as a function of
     /// RoiCuts::energies()
     std::vector<double> m_exposure;
     
     /// Energies for the diffuse exposure calculation
     std::vector<double> m_energies;
     
     
     /// FIXME, what exactly are these?
     void getExposureArrays(double emin, double emax,
			    std::vector<double> & energies,
			    std::vector<double> & exposures,
			    size_t nee = 0) const;
     
     /// FIXME, what exactly are these?
     void getExposureSubArrays(double emin, double emax,
			       std::vector<double> & energies,
			       std::vector<double> & exposures) const;
     
     /// FIXME, what exactly are these?
     void getExposureValues(const std::vector<double> & energies,
			    std::vector<double> & exposures) const;
     
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
	 : optimizers::Function("FuncDeriv", func.getNumParams(), ""),
	   m_func(func), m_parName(parName) {}
     protected:
       virtual double value(const optimizers::Arg & x) const {
	 return m_func.derivByParam(x, m_parName);
       }
       virtual double derivByParamImp(const optimizers::Arg &,
				      const std::string &) const {
	 throw std::runtime_error("FluxDeriv::deriveByParam not implemented");
       }
       
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
       EnergyFlux(const optimizers::Function & func)  
         : optimizers::Function("EnergyFlux", func.getNumParams(), ""),
           m_func(func) {}
     protected:
       virtual double value(const optimizers::Arg & x) const {
         double energy(dynamic_cast<const optimizers::dArg &>(x).getValue());
         return energy*m_func(x);
       }
       virtual double derivByParamImp(const optimizers::Arg & x,
				      const std::string & parname) const {
         double energy(dynamic_cast<const optimizers::dArg &>(x).getValue());
         return energy*m_func.derivByParam(x, parname);
       }
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
         : optimizers::Function("EnergyFluxDeriv", func.getNumParams(), ""),
           m_func(func), m_parName(parName) {}
     protected:
       virtual double value(const optimizers::Arg & x) const {
         double energy(dynamic_cast<const optimizers::dArg &>(x).getValue());
         return energy*m_func.derivByParam(x, m_parName);
       }
       virtual double derivByParamImp(const optimizers::Arg &,
				      const std::string &) const {
         throw std::runtime_error("EnergyFluxDeriv::derivByParam: "
                                  "not implemented");
       }
       virtual Function * clone() const {
         return 0;
       }
     private:
       const optimizers::Function & m_func;
       std::string m_parName;
     };
     

     /** 
      * Template function to compute integral over energy.

      Typically this takes one of the three classes above as input 
      **/
     template<typename Functor>
     double computeEnergyIntegral(const Functor & func, 
				  double emin, double emax, size_t npts) const {
       std::vector<double> energies;
       energies.reserve(npts);
       double estep(std::log(emax/emin)/float(npts-1));
       for (size_t k(0); k < npts; k++) {
         energies.push_back(emin*std::exp(estep*k));
       }
       return computeEnergyIntegral(func, energies);
     }

     
     /** 
      * Template function to compute integral over energy.

      Typically this takes one of the three classes above as input 
      **/
     template<typename Functor>
     double computeEnergyIntegral(const Functor & func, 
				  const std::vector<double> & energies) const {
       std::vector<double> integrand;
       integrand.reserve(energies.size());
       for (size_t k(0); k < energies.size(); k++) {
         optimizers::dArg arg(energies.at(k));
         integrand.push_back(func(arg)*angularIntegral(energies.at(k)));
       }
       bool useLog;
       TrapQuad fluxIntegral(energies, integrand, useLog=true);
       return fluxIntegral.integral();
     }

   };

} // namespace Likelihood

#endif // Likelihood_Source_h
