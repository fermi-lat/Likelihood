/** 
 * @file Source.h
 * @brief Source base class declaration
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Source.h,v 1.21 2003/08/13 18:01:15 jchiang Exp $
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
 * @brief Abstract base class for gamma-ray sources.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Source.h,v 1.21 2003/08/13 18:01:15 jchiang Exp $
 */

class Source {

public:
    
   Source() {m_name = "";}
   Source(const Source &rhs);
   virtual ~Source() {}

   //! returns photons/cm^2-s-sr-MeV having been convolved through
   //! the LAT instrument response
   virtual double fluxDensity(const Event &evt) const = 0;

   //! derivatives of fluxDensity wrt model Parameters
   virtual double fluxDensityDeriv(const Event &evt, 
                                   const std::string &paramName) const = 0;

   //! predicted number of photons given RoiCuts and ScData
   virtual double Npred() = 0;

   //! derivative of Npred wrt named Parameter
   virtual double NpredDeriv(const std::string &paramName) = 0;

   //! access unique source identifier
   void setName(const std::string &name) {m_name = name;}
   std::string getName() const {return m_name;}

   //! return a reference to the m_functions map (NB: not const!)
   typedef std::map<std::string, optimizers::Function *> FuncMap;
   FuncMap & getSrcFuncs() {return m_functions;}

   virtual void setDir(double ra, double dec, bool updateExposure = true) = 0;
   virtual void setDir(const astro::SkyDir &dir, 
                       bool updateExposure = true) = 0;

   virtual void setSpectrum(optimizers::Function *) = 0;
                       
   //! clone function, with default
   virtual Source *clone() const {return 0;}

   //! return the Source type (e.g., Diffuse vs Point)
   std::string getType() {return m_srcType;}

protected:

   //! Source name
   std::string m_name;

   //! Source type
   std::string m_srcType;

   //! map of Functions describing this source
   FuncMap m_functions;

};

} // namespace Likelihood

#endif // Likelihood_Source_h
