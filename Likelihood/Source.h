/** 
 * @file Source.h
 * @brief Source base class declaration
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Source.h,v 1.30 2005/02/28 18:38:46 jchiang Exp $
 */

#ifndef Likelihood_Source_h
#define Likelihood_Source_h

#include <iostream>
#include <map>
#include "astro/SkyDir.h"
#include "optimizers/Function.h"
#include "Likelihood/Event.h"

namespace Likelihood {

/** 
 * @class Source
 *
 * @brief Base class for gamma-ray sources.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Source.h,v 1.30 2005/02/28 18:38:46 jchiang Exp $
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

   /// derivatives of fluxDensity wrt model Parameters
   virtual double fluxDensityDeriv(const Event &evt, 
                                   const std::string &paramName) const = 0;

   virtual double fluxDensityDeriv(double inclination, double phi, 
                                   double energy, 
                                   const astro::SkyDir & appDir,
                                   int evtType, 
                                   const std::string & paramName) const = 0;

   /// predicted number of photons given RoiCuts and ScData
   virtual double Npred() = 0;

   /// derivative of Npred wrt named Parameter
   virtual double NpredDeriv(const std::string &paramName) = 0;

   /// @return predicted number of counts within a specified energy range
   virtual double Npred(double emin, double emax) = 0;

   /// access unique source identifier
   void setName(const std::string & name) {
      m_name = name;
   }

   const std::string & getName() const {
      return m_name;
   }

   /// @return a reference to the m_functions map (NB: not const!)
   typedef std::map<std::string, optimizers::Function *> FuncMap;
   FuncMap & getSrcFuncs() {
      return m_functions;
   }

   virtual void setSpectrum(optimizers::Function *) = 0;
                       
   virtual Source * clone() const {
      return 0;
   }

   /// @return the Source type (e.g., Diffuse vs Point)
   const std::string & getType() const {
      return m_srcType;
   }

   /// Integrate the product of the source spectrum with the given
   /// SourceMap pixel values.
   virtual double pixelCounts(double emin, double emax, 
                              double wtMin, double wtMax) const = 0;

   virtual double pixelCountsDeriv(double emin, double emax, 
                                   double wtMin, double wtMax,
                                   const std::string & paramName) const = 0;

protected:

   /// Source name
   std::string m_name;

   /// Source type
   std::string m_srcType;

   /// map of Functions describing this source
   FuncMap m_functions;

   /// flag to indicate if energy dispersion is to be used.
   bool m_useEdisp;

};

} // namespace Likelihood

#endif // Likelihood_Source_h
