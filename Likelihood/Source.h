/** 
 * @file Source.h
 * @brief Source base class declaration
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Source.h,v 1.33 2005/03/02 04:51:10 jchiang Exp $
 */

#ifndef Likelihood_Source_h
#define Likelihood_Source_h

#include <iostream>
#include <map>

#include "optimizers/Function.h"

#include "Likelihood/Event.h"

namespace astro {
   class SkyDir;
}

namespace Likelihood {

/** 
 * @class Source
 *
 * @brief Base class for gamma-ray sources.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Source.h,v 1.33 2005/03/02 04:51:10 jchiang Exp $
 */

class Source {

public:
    
   Source();

   Source(const Source &rhs);

   virtual ~Source() {}

   /// @return photons/cm^2-s-sr-MeV having been convolved through
   /// the LAT instrument response
   virtual double fluxDensity(const Event &evt) const = 0;

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
                              int evtType) const = 0;

   /// Derivatives of fluxDensity wrt model Parameters
   virtual double fluxDensityDeriv(const Event &evt, 
                                   const std::string &paramName) const = 0;

   virtual double fluxDensityDeriv(double inclination, double phi, 
                                   double energy, const astro::SkyDir & appDir,
                                   int evtType, 
                                   const std::string & paramName) const = 0;

   /// Predicted number of photons.
   virtual double Npred() = 0;

   /// Derivative of Npred wrt named Parameter
   virtual double NpredDeriv(const std::string &paramName) = 0;

   /// Predicted number of counts within a specified energy range.
   virtual double Npred(double emin, double emax) = 0;

   /// Set the Function used for modeling the source spectrum.
   virtual void setSpectrum(optimizers::Function *) = 0;
                       
   virtual void setName(const std::string & name) {
      m_name = name;
   }

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
                              double wtMin, double wtMax) const = 0;

   virtual double pixelCountsDeriv(double emin, double emax, 
                                   double wtMin, double wtMax,
                                   const std::string & paramName) const = 0;

   virtual const std::vector<double> & exposure() const = 0;

protected:

   /// A unique source name.
   std::string m_name;

   /// Source type --- "Diffuse" or "Point".
   std::string m_srcType;

   /// Flag to indicate if energy dispersion is to be used.
   bool m_useEdisp;

   /// Map of Functions describing this source.
   std::map<std::string, optimizers::Function *> m_functions;

};

} // namespace Likelihood

#endif // Likelihood_Source_h
