#ifndef Source_h
#define Source_h

#include <map>
#include "astro/SkyDir.h"
#include "../Likelihood/Function.h"
#include "../Likelihood/Event.h"

namespace Likelihood {

/** 
 * @class Source
 *
 * @brief Abstract base class for gamma-ray sources.
 *
 * @author J. Chiang
 *    
 * $Header: 
 */

class Source {

public:
    
   Source(){m_name = "";};
   Source(const Source &rhs) {m_name = rhs.m_name;};
   virtual ~Source(){};

   //! returns photons/cm^2-s-sr-GeV having been convolved through
   //! the LAT instrument response
   virtual double fluxDensity(double energy, double time,
                              const astro::SkyDir &dir) const = 0;

   double fluxDensity(const Event &evt) const
      {return fluxDensity(evt.getEnergy(), evt.getArrTime(), evt.getDir());}

   //! derivatives of fluxDensity wrt model Parameters
   virtual double fluxDensityDeriv(double energy, double time,
                                   const astro::SkyDir &dir,
                                   std::string &paramName) const = 0;

   double fluxDensityDeriv(const Event &evt, std::string &paramName) const
      {return fluxDensityDeriv(evt.getEnergy(), evt.getArrTime(), 
                               evt.getDir(), paramName);}

   //! predicted number of photons given RoiCuts and ScData
   virtual double Npred() = 0;

   //! derivative of Npred wrt named Parameter
   virtual double NpredDeriv(const std::string &paramName) = 0;

   //! access unique source identifier
   void setName(const std::string &name) {m_name = name;};
   std::string getName() const {return m_name;};

   //! return a reference to the m_functions map (NB: not const!)
   typedef std::map<std::string, Function *> FuncMap;
   FuncMap & getSrcFuncs() {return m_functions;}

protected:

   //! source name
   std::string m_name;

   //! map of Functions describing this source
   FuncMap m_functions;

};

} // namespace Likelihood

#endif // Source_h
