/** 
 * @file ResponseFunctions.h
 * @brief A singleton class to contain the instrument response functions.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ResponseFunctions.h,v 1.13 2005/02/27 06:42:24 jchiang Exp $
 */

#ifndef Likelihood_ResponseFunctions_h
#define Likelihood_ResponseFunctions_h

#include <map>

#include "irfInterface/Irfs.h"

namespace astro {
   class SkyDir;
}

namespace Likelihood {

/** 
 * @class ResponseFunctions
 *
 * @brief This class provides global access to a map of pointers to
 * irfInterface::Irfs objects.  These pointers are indexed by event
 * type, given as an integer; a map is used since the indices need not
 * be contiguous.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ResponseFunctions.h,v 1.13 2005/02/27 06:42:24 jchiang Exp $
 */

class ResponseFunctions {
    
public:
    
   virtual ~ResponseFunctions() {}

   static ResponseFunctions * instance();

   /// Return the total instrument response 
   /// (= effective area*PSF*energy dispersion).
   /// @param energy True photon energy (MeV).
   /// @param appEnergy Measured energy (MeV).
   /// @param zAxis Direction of spacecraft z-axis.
   /// @param xAxis Direction of spacecraft x-axis.
   /// @param srcDir Assumed source (i.e., true photon) direction.
   /// @param appDir Apparent photon direction.
   /// @param type Event type identifying which set of IRFs to use.
   ///        Presently, 0 = Front part of LAT, 1 = Back part of LAT.
   ///        2 = Combined (GLAST25 only). 
   ///        (@todo These IDs need to be rationalized and coordinated 
   ///        with the irfInterface package.)
   static double totalResponse(double energy, double appEnergy,
                               const astro::SkyDir & zAxis,
                               const astro::SkyDir & xAxis,
                               const astro::SkyDir & srcDir,
                               const astro::SkyDir & appDir,
                               int type);

   static double totalResponse(double inclination, double phi, 
                               double energy, double appEnergy, 
                               double separation, int evtType);

   static void setRespPtrs(std::map<unsigned int, irfInterface::Irfs *> 
                           &respPtrs) {s_respPtrs = respPtrs;}

   static void addRespPtr(unsigned int key,
                          irfInterface::Irfs *respPtr) {
      s_respPtrs[key] = respPtr;
   }

   static void deleteRespPtr(unsigned int key) {
      if (s_respPtrs.count(key)) {
         delete s_respPtrs[key];
         s_respPtrs[key] = 0;
      }
   }

   irfInterface::Irfs * respPtr(unsigned int eventType);

   std::map<unsigned int, irfInterface::Irfs *>::const_iterator begin() const {
      return s_respPtrs.begin();
   }

   std::map<unsigned int, irfInterface::Irfs *>::const_iterator end() const {
      return s_respPtrs.end();
   }

   /// Whether or not energy dispersion is to be considered.
   static const bool & useEdisp() {return s_useEdisp;}

   static void setEdispFlag(bool useEdisp) {s_useEdisp = useEdisp;}

   static const std::string & respName() {return s_respName;}

   static void setRespName(const std::string & respName) {
      s_respName = respName;
   }

protected:

   ResponseFunctions() {}

private:

   static ResponseFunctions * s_instance;

   static std::map<unsigned int, irfInterface::Irfs *> s_respPtrs;

   static bool s_useEdisp;

   static std::string s_respName;

};

} // namespace Likelihood

#endif // Likelihood_ResponseFunctions_h
