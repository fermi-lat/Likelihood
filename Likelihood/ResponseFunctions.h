/** 
 * @file ResponseFunctions.h
 * @brief A class to contain the instrument response functions.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ResponseFunctions.h,v 1.14 2005/03/02 22:55:00 jchiang Exp $
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
 * @brief This class provides access to a map of pointers to
 * irfInterface::Irfs objects.  These pointers are indexed by event
 * type, given as an integer; a map is used since the indices need not
 * be contiguous.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ResponseFunctions.h,v 1.14 2005/03/02 22:55:00 jchiang Exp $
 */

class ResponseFunctions {
    
public:
    
   ResponseFunctions() : m_useEdisp(false), m_respName("") {}

   virtual ~ResponseFunctions() {}

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
   double totalResponse(double energy, double appEnergy,
                        const astro::SkyDir & zAxis,
                        const astro::SkyDir & xAxis,
                        const astro::SkyDir & srcDir,
                        const astro::SkyDir & appDir,
                        int type) const;

   double totalResponse(double inclination, double phi, 
                        double energy, double appEnergy, 
                        double separation, int evtType) const;
   
   void setRespPtrs(std::map<unsigned int, irfInterface::Irfs *> 
                    &respPtrs) {
      m_respPtrs = respPtrs;
   }

   void addRespPtr(unsigned int key,
                   irfInterface::Irfs *respPtr) {
      m_respPtrs[key] = respPtr;
   }

   void deleteRespPtr(unsigned int key) {
      if (m_respPtrs.count(key)) {
         delete m_respPtrs[key];
         m_respPtrs[key] = 0;
      }
   }

   irfInterface::Irfs * respPtr(unsigned int eventType) const;

   std::map<unsigned int, irfInterface::Irfs *>::const_iterator begin() const {
      return m_respPtrs.begin();
   }

   std::map<unsigned int, irfInterface::Irfs *>::const_iterator end() const {
      return m_respPtrs.end();
   }

   /// Whether or not energy dispersion is to be considered.
   const bool & useEdisp() const {
      return m_useEdisp;
   }

   void setEdispFlag(bool useEdisp) {
      m_useEdisp = useEdisp;
   }

   const std::string & respName() const {
      return m_respName;
   }

   void setRespName(const std::string & respName) {
      m_respName = respName;
   }

private:

   std::map<unsigned int, irfInterface::Irfs *> m_respPtrs;

   bool m_useEdisp;

   std::string m_respName;

};

} // namespace Likelihood

#endif // Likelihood_ResponseFunctions_h
