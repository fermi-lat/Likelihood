/** 
 * @file ResponseFunctions.h
 * @brief A singleton class to contain the instrument response functions.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ResponseFunctions.h,v 1.4 2004/02/23 22:16:33 jchiang Exp $
 */

#ifndef Likelihood_ResponseFunctions_h
#define Likelihood_ResponseFunctions_h

#include <map>

#include "latResponse/Irfs.h"

namespace astro {
   class SkyDir;
}

namespace Likelihood {

/** 
 * @class ResponseFunctions
 *
 * @brief This class provides global access to a map of pointers to
 * latResponse::Irfs objects.  These pointers are indexed by event
 * type, given as an integer; a map is used since the indices need not
 * be contiguous.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ResponseFunctions.h,v 1.4 2004/02/23 22:16:33 jchiang Exp $
 */

class ResponseFunctions {
    
public:
    
   virtual ~ResponseFunctions() {}

   static ResponseFunctions * instance();

   /// Return the total instrument response 
   /// (= effective area*PSF*energy dispersion).
   /// @param time MET (seconds).  These are the same time units used
   ///        by the spacecraft data class ScData.
   /// @param energy True photon energy (MeV).
   /// @param appEnergy Measured energy (MeV).
   /// @param srcDir Assumed source (i.e., true photon) direction.
   /// @param appDir Apparent photon direction.
   /// @param type Event type identifying which set of IRFs to use.
   ///        Presently, 0 = Front part of LAT, 1 = Back part of LAT.
   ///        2 = Combined (GLAST25 only). 
   ///        (@todo These IDs need to be rationalized and coordinated 
   ///        with the latResponse package.)
   static double totalResponse(double time,
                               double energy, double appEnergy,
                               const astro::SkyDir &srcDir,
                               const astro::SkyDir &appDir,
                               int type);

   static void setRespPtrs(std::map<unsigned int, latResponse::Irfs *> 
                           &respPtrs) {s_respPtrs = respPtrs;}

   static void addRespPtr(unsigned int key,
                          latResponse::Irfs *respPtr) 
      {s_respPtrs[key] = respPtr;}

   static void deleteRespPtr(unsigned int key) {
      if (s_respPtrs.count(key)) {
         delete s_respPtrs[key];
         s_respPtrs[key] = 0;
      }
   }

   latResponse::Irfs * respPtr(unsigned int eventType);

   std::map<unsigned int, latResponse::Irfs *>::iterator begin()
      {return s_respPtrs.begin();}

   std::map<unsigned int, latResponse::Irfs *>::iterator end()
      {return s_respPtrs.end();}


protected:

   ResponseFunctions() {}

private:

   static ResponseFunctions * s_instance;

   static std::map<unsigned int, latResponse::Irfs *> s_respPtrs;

};

} // namespace Likelihood

#endif // Likelihood_ResponseFunctions_h
