/** 
 * @file RoiCuts.h
 * @brief Declaration for RoiCuts class
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/RoiCuts.h,v 1.26 2005/03/01 01:06:52 jchiang Exp $
 */

#ifndef Likelihood_RoiCuts_h
#define Likelihood_RoiCuts_h

#include <cmath>

#include <string>
#include <utility>
#include <vector>

#include "astro/SkyDir.h"

#include "irfInterface/AcceptanceCone.h"

#include "dataSubselector/Cuts.h"
#include "dataSubselector/GtiCut.h"
#include "dataSubselector/RangeCut.h"
#include "dataSubselector/SkyConeCut.h"

namespace tip {
   class Header;
}

namespace Likelihood {

class Event;

/** 
 * @class RoiCuts
 *
 * @brief NTuple Singleton class to represent Region-of-interest cuts.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/RoiCuts.h,v 1.26 2005/03/01 01:06:52 jchiang Exp $
 */

class RoiCuts {

public:

   static RoiCuts * instance();

   /// Access to the accepted time intervals.  An event time is valid
   /// if it lies within the union of these intervals.
   void getTimeCuts(std::vector< std::pair<double, double> > &timeCuts) const {
      timeCuts = s_timeCuts;
   }

   std::pair<double, double> getEnergyCuts() const {
      return std::make_pair(s_eMin, s_eMax);
   }

   const irfInterface::AcceptanceCone & extractionRegion() const {
      return s_roiCone;
   }

   void getRaDec(double & ra, double & dec) const {
      ra = s_roiCone.center().ra();
      dec = s_roiCone.center().dec();
   }

   double getMuZenMax() const {
      return s_muZenMax;
   }

   /// Set all cuts (includes reset of time cuts)
   void setCuts(double ra = 193.98, double dec = -5.82, 
                double roi_radius = 20.,
                double emin = 30., double emax = 3.1623e5,
                double tmin = 0., double tmax = 1e12,
                double muZenMax = -1.);

   /// Read from the DSS keywords in the eventFile. (Series of event
   /// files are required to have the same DSS keywords.)
   void readCuts(const std::string & eventFile,
                 const std::string & ext="EVENTS",
                 bool strict=true);

   /// Apply these cuts to an Event
   bool accept(const Event &) const;

   /// Write DSS keywords to a FITS header
   void writeDssKeywords(tip::Header & header) const;

   /// Write the GTI extension to the FITS file
   void writeGtiExtension(const std::string & filename) const;

   /// A logrithmically spaced vector of energies from the minimum
   /// energy to the maximum energy.
   const std::vector<double> & energies() const {
      return m_energies;
   }

protected:

   RoiCuts() : m_energyCut(0), m_skyConeCut(0) {
      makeEnergyVector();
   }

   ~RoiCuts() {
      for (int i = m_gtiCuts.size()-1; i > -1; i--) {
         delete m_gtiCuts.at(i);
      }
      for (int i = m_timeCuts.size()-1; i > -1; i--) {
         delete m_timeCuts.at(i);
      }
      delete m_skyConeCut;
      delete m_energyCut;
   }

private:

   static RoiCuts * s_instance;

   static dataSubselector::Cuts * s_cuts;

   /// cuts on photon "MET" arrival times in seconds; 
   /// this vector of pairs specify time intervals for event acceptance;
   /// the *intersection* of these intervals will be used
   typedef std::pair<double, double> timeInterval; // this will be generalized
   static std::vector<timeInterval> s_timeCuts;

   /// minimum and maximum energies in MeV,
   static double s_eMin;
   static double s_eMax;

   /// A vector of logrithmically-spaced energies between s_eMin and s_eMax.
   std::vector<double> m_energies;

   /// The acceptance cone or sky extraction region.
   static irfInterface::AcceptanceCone s_roiCone;

   /// cosine of the maximum Zenith angle
   static double s_muZenMax;

   dataSubselector::RangeCut * m_energyCut;
   dataSubselector::SkyConeCut * m_skyConeCut;
   std::vector<dataSubselector::RangeCut *> m_timeCuts;
   std::vector<dataSubselector::GtiCut *> m_gtiCuts;

   /// Set additional time cuts
   void addTimeInterval(double tmin, double tmax);

   /// Create the m_energies vector.
   void makeEnergyVector(int nee=100);

   void sortCuts(bool strict=true);
   void setRoiData();
};

} // namespace Likelihood

#endif // Likelihood_RoiCuts_h

