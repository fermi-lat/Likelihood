#ifndef Source_h
#define Source_h

#include "astro/SkyDir.h"

namespace Likelihood {

/** 
 * @class Source
 *
 * @brief Base class for gamma-ray sources.
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
   virtual double fluxDensity(const double energy, const double time,
			      const astro::SkyDir &dir) const = 0;

   //! access unique source identifier
   void setName(const std::string &name) {m_name = name;};
   std::string getName() const {return m_name;};

private:

   //! source name
   std::string m_name;

};

} // namespace Likelihood

#endif // Source_h
