#ifndef Source_h
#define Source_h

#include <map>
#include "astro/SkyDir.h"
#include "../Likelihood/Function.h"
#include "../Likelihood/ScData.h"

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
    
   Source() {
      ScData *scData = ScData::instance();
      m_name = "";}
   Source(const Source &rhs) {m_name = rhs.m_name;};
   virtual ~Source(){};

   //! returns photons/cm^2-s-sr-GeV having been convolved through
   //! the LAT instrument response
   virtual double fluxDensity(const double energy, const double time,
			      const astro::SkyDir &dir) const = 0;

   //! the predicted number of photons for a given ROI
   //! the vector<double> is a place holder for an ROI_cuts class object
   virtual double Npred(const std::vector<double>) const = 0;

   //! access unique source identifier
   void setName(const std::string &name) {m_name = name;};
   std::string getName() const {return m_name;};

   //! return a reference to the m_functions map (NB: not const!)
   std::map<std::string, Function *> & getSrcFuncs() {return m_functions;}

protected:

   //! spacecraft data for all Sources to share
   ScData * scData;

   //! source name
   std::string m_name;

   //! map of Functions describing this source
   typedef std::map<std::string, Function *> FuncMap;
   FuncMap m_functions;

};

} // namespace Likelihood

#endif // Source_h
