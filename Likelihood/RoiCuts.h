/** 
 * @file RoiCuts.h
 * @brief Declaration for RoiCuts class
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/RoiCuts.h,v 1.14 2004/02/20 00:02:03 jchiang Exp $
 */

#ifndef Likelihood_RoiCuts_h
#define Likelihood_RoiCuts_h

#include <vector>
#include <string>
#include <utility>
#include <cmath>

#include "xml/Dom.h"

#include "astro/SkyDir.h"

#include "latResponse/AcceptanceCone.h"

namespace Likelihood {

class Event;

/** 
 * @class RoiCuts
 *
 * @brief NTuple Singleton class to represent Region-of-interest cuts.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/RoiCuts.h,v 1.14 2004/02/20 00:02:03 jchiang Exp $
 */

class RoiCuts {

public:

   ~RoiCuts(){}

   static RoiCuts * instance();

   /// Access to the cuts
   void getTimeCuts(std::vector< std::pair<double, double> > &tLimVec) const
      {tLimVec = s_tLimVec;}

   std::pair<double, double> getEnergyCuts() const
      {return std::make_pair(s_eMin, s_eMax);}

   const latResponse::AcceptanceCone &extractionRegion() const
      {return s_roiCone;}

   static void getRaDec(double &ra, double &dec) {
      ra = s_roiCone.center().ra();
      dec = s_roiCone.center().dec();
   }

   double getMuZenMax() {return s_muZenMax;}

   /// Methods to allow cuts to be specified

   /// Set additional time cuts
   static void addTimeInterval(double tmin, double tmax);

   /// Set all cuts (includes reset of time cuts)
   static void setCuts(double ra = 193.98, double dec = -5.82, 
                       double roi_radius = 20.,
                       double emin = 30., double emax = 3.1623e5,
                       double tmin = 0., double tmax = 1e12,
                       double muZenMax = -1.);

   /// Read from xml file
   static void setCuts(std::string xmlFile);

   /// Write to xml file
   void writeXml(std::string xmlFile, const std::string &roiTitle="");

   /// Write to a stream.
   void writeXml(std::ostream & ostr, const std::string &roiTitle="",
                 bool pretty=false);

   /// Apply these cuts to an Event
   bool accept(const Event &);

protected:

   RoiCuts(){}

private:

   static RoiCuts * s_instance;

   /// cuts on photon "MET" arrival times in seconds; 
   /// this vector of pairs specify time intervals for event acceptance;
   /// the *intersection* of these intervals will be used
   typedef std::pair<double, double> timeInterval; // this will be generalized
   static std::vector<timeInterval> s_tLimVec;

   /// minimum and maximum energies in MeV,
   static double s_eMin;
   static double s_eMax;

   /// The acceptance cone or sky extraction region.
   static latResponse::AcceptanceCone s_roiCone;

   /// cosine of the maximum Zenith angle
   static double s_muZenMax;

   /// Create a root DomElement for the current set of cuts.
   DomElement * rootDomElement(const std::string & roiTitle);

};

} // namespace Likelihood
#endif // Likelihood_RoiCuts_h
