/** 
 * @file RoiCuts.h
 * @brief Declaration for RoiCuts class
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/RoiCuts.h,v 1.38 2013/08/26 22:55:31 jchiang Exp $
 */

#ifndef Likelihood_RoiCuts_h
#define Likelihood_RoiCuts_h

#include <cmath>

#include <string>
#include <utility>
#include <vector>

#include "astro/SkyDir.h"

#include "irfInterface/AcceptanceCone.h"

#include "dataSubselector/BitMaskCut.h"
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
 * @brief Class to represent Region-of-interest cuts.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/RoiCuts.h,v 1.38 2013/08/26 22:55:31 jchiang Exp $
 */

class RoiCuts {

public:

   RoiCuts() : m_cuts(0), m_eMin(20.), m_eMax(2e5), 
               m_energyCut(0), m_skyConeCut(0), m_minTime(0), m_maxTime(0) {
      makeEnergyVector();
   }

   ~RoiCuts() {
      delete m_cuts;
   }

   /// Access to the time range cut intervals.  An event time is valid
   /// only if it lies within the union of these intervals.
   void getTimeCuts(std::vector< std::pair<double, double> > &timeCuts) const;

   void getGtis(std::vector< std::pair<double, double> > & gtis) const;

   std::pair<double, double> getEnergyCuts() const {
      return std::make_pair(m_eMin, m_eMax);
   }

#ifndef SWIG
   const irfInterface::AcceptanceCone & extractionRegion() const {
      return m_roiCone;
   }
#endif //SWIG

   std::vector<double> roiCone() const {
      std::vector<double> my_vec;
      my_vec.push_back(m_roiCone.center().ra());
      my_vec.push_back(m_roiCone.center().dec());
      my_vec.push_back(m_roiCone.radius());
      return my_vec;
   }

   void getRaDec(double & ra, double & dec) const {
      ra = m_roiCone.center().ra();
      dec = m_roiCone.center().dec();
   }

   double getMuZenMax() const {
      return m_muZenMax;
   }

   /// Set all cuts (includes reset of time cuts)
   void setCuts(double ra = 193.98, double dec = -5.82, 
                double roi_radius = 20.,
                double emin = 30., double emax = 3.1623e5,
                double tmin = 0., double tmax = 1e12,
                double muZenMax = -1., bool reset_tlims=false);

   /// Read from the DSS keywords in a single eventFile.
   void readCuts(const std::string & eventFile,
                 const std::string & ext="EVENTS",
                 bool strict=true);

   /// Read from the DSS keywords from several eventFiles. These
   /// files are required to have the same DSS keywords, but the
   /// GTIs need not be the same. The GTIs from the various files 
   /// will be concatenated into a single Gti object.
   void readCuts(const std::vector<std::string> & eventFiles,
                 const std::string & ext="EVENTS",
                 bool strict=true);

   /// Apply these cuts to an Event
   bool accept(const Event &) const;

   /// Write DSS keywords to a FITS header
   void writeDssKeywords(tip::Header & header) const;

   /// Write DSS time-related keywords to a FITS header
   void writeDssTimeKeywords(tip::Header & header) const;

   /// Write the GTI extension to the FITS file
   void writeGtiExtension(const std::string & filename) const;

   /// A logrithmically spaced vector of energies from the minimum
   /// energy to the maximum energy.
   const std::vector<double> & energies() const {
      return m_energies;
   }

   double minTime() const {
      return m_minTime;
   }

   double maxTime() const {
      return m_maxTime;
   }

   void setIrfsVersion(const std::string & irfsName);

   void setBitMaskCut(dataSubselector::BitMaskCut * candidateCut) {
      m_cuts->setBitMaskCut(candidateCut);
   }

private:

   dataSubselector::Cuts * m_cuts;

   /// minimum and maximum energies in MeV,
   double m_eMin;
   double m_eMax;

   /// A vector of logrithmically-spaced energies between m_eMin and m_eMax.
   std::vector<double> m_energies;

   /// The acceptance cone or sky extraction region.
   irfInterface::AcceptanceCone m_roiCone;

   /// cosine of the maximum Zenith angle
   double m_muZenMax;

   dataSubselector::RangeCut * m_energyCut;
   dataSubselector::SkyConeCut * m_skyConeCut;
   std::vector<dataSubselector::RangeCut *> m_timeRangeCuts;
   std::vector<dataSubselector::GtiCut *> m_gtiCuts;

   /// Minimum and maximum spacecraft times to be considered based
   /// on time range cuts and GTIs
   double m_minTime;
   double m_maxTime;

   /// Add a time range cut.
   void addTimeInterval(double tmin, double tmax);

   /// Create the m_energies vector.
   void makeEnergyVector(int nee=100);

   void sortCuts(bool strict=true);
   void setRoiData();
};

} // namespace Likelihood

#endif // Likelihood_RoiCuts_h

