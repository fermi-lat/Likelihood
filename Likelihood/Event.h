/** @file Event.h
 * @brief Event class declaration
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef Event_h
#define Event_h

#include <vector>
#include <string>
#include <map>
#include <utility>

#include "astro/SkyDir.h"

namespace Likelihood {

/** 
 * @class Event
 *
 * @brief A gamma-ray event --- apparent direction, energy, arrival time, etc.
 *
 * @author J. Chiang
 *    
 * $Header$
 */

class Event {
    
public:

   Event(){};
   Event(double ra, double dec, double energy, double time, 
         double sc_ra, double sc_dec, double muZenith);
   Event(const Event &);
   virtual ~Event(){}

   astro::SkyDir getDir() const {return m_appDir;}
   astro::SkyDir getScDir() const {return m_scDir;}

   //! some useful accessor functions
   double getEnergy() const {return m_energy;}
   double getArrTime() const {return m_arrTime;}
   double getMuZenith() const {return m_muZenith;}

   //! separation in units of radians
   double getSeparation(const astro::SkyDir &dir) const 
      {return m_appDir.SkyDir::difference(dir);}
    
private:

   //! apparent direction, energy, arrival time, and cosine(zenith angle)
   astro::SkyDir m_appDir;
   double m_energy;
   double m_arrTime;
   double m_muZenith;
   
   //! spacecraft info at event arrival time
   astro::SkyDir m_scDir;

   //! response function data, unique to each event, and 
   //! comprising an energy redistribution function...
   typedef std::vector< std::pair<double, double> > diffuse_response;

   //! for uniform diffuse extragalactic emission,
   diffuse_response m_respEg;

   //! for the diffuse Galactic model,
   diffuse_response m_respGal;

   //! and for any number of diffuse sources,
   std::map<std::string, diffuse_response> m_respDiffuseSrcs;

   //! something to compute the event response function data
   void computeResponse()
      {std::cout << "Computing event response..." << std::endl;}
        
};

} // namespace Likelihood
#endif // Event_h
