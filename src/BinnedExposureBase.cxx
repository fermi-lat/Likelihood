/**
 * @file BinnedExposureBase.cxx
 * @brief Integral of effective area over time for the entire sky at
 * various energies.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/echarles/healpix_changes/Likelihood/src/BinnedExposureBase.cxx,v 1.2 2015/03/03 05:59:59 echarles Exp $
 */

#include <cmath>
#include <cstdio>

#include <algorithm>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "st_stream/StreamFormatter.h"

#include "st_app/AppParGroup.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Image.h"

#include "astro/SkyProj.h"

#include "st_facilities/Util.h"

#include "Likelihood/BinnedExposureBase.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/Observation.h"


namespace Likelihood {

std::vector<double>::const_iterator
BinnedExposureBase::findNearest(const std::vector<double> & xx, double x, double tol) {
  std::vector<double>::const_iterator ix = std::find(xx.begin(),
						     xx.end(), x);
  if (ix == xx.end()) { // no exact match, so look for nearest                                                       
    for (ix = xx.begin(); ix != xx.end(); ++ix) {
      if (fracDiff(x, *ix) < tol) {
	return ix;
      }
    }
    std::ostringstream what;
    what << "BinnedExposureBase::operator(): The energy " << x
	 << " is not available.\nHere are the relevant energies:\n";
    for (size_t i(0); i < xx.size(); i++) {
      what << xx.at(i) << "\n";
    }
    throw std::runtime_error(what.str());
  }
  return ix;  // return the exact match                                                                              
}

BinnedExposureBase::BinnedExposureBase() : m_observation(0), m_proj(0), 
					   m_costhmin(-1), m_costhmax(1),
					   m_enforce_boundaries(false) {}

BinnedExposureBase::BinnedExposureBase(const Observation & observation,
				       bool useEbounds,
				       const st_app::AppParGroup * pars)
   : m_observation(&observation), m_proj(0), m_costhmin(-1), m_costhmax(1),
     m_enforce_boundaries(false),
     m_allSky(false){
   if (pars) {
      setCosThetaBounds(*pars);
   }
}

BinnedExposureBase::BinnedExposureBase(const std::vector<double> & energies,
				       const Observation & observation,
				       const st_app::AppParGroup * pars) 
   : m_energies(energies), m_observation(&observation), m_proj(0),
     m_costhmin(-1), m_costhmax(1),m_enforce_boundaries(false),m_allSky(false)  {
   if (pars) {
      setCosThetaBounds(*pars);
   } 
}

BinnedExposureBase::BinnedExposureBase(const std::string & filename) 
   : m_observation(0), m_proj(0), m_costhmin(-1), m_costhmax(1),
     m_enforce_boundaries(false),m_allSky(false) {

   std::auto_ptr<const tip::Table>
    energies(tip::IFileSvc::instance().readTable(filename, "Energies"));

   m_energies.clear();
   tip::Table::ConstIterator it = energies->begin();
   tip::ConstTableRecord & row = *it;
   for ( ; it != energies->end(); ++it) {
      double value;
      row["Energy"].get(value);
      m_energies.push_back(value);
   }
}

BinnedExposureBase::~BinnedExposureBase() {
  delete m_proj;
}


void BinnedExposureBase::setCosThetaBounds(const st_app::AppParGroup & pars) {
   double thmin = pars["thmin"];
   if (thmin > 0) {
      m_costhmax = std::cos(thmin*M_PI/180.);
   }
   double thmax = pars["thmax"];
   if (thmax < 180.) {
      m_costhmin = std::cos(thmax*M_PI/180.);
   }
}

double BinnedExposureBase::Aeff::value(double cosTheta, double phi) const {
   if (cosTheta < m_costhmin || cosTheta > m_costhmax) {
      return 0;
   }
   return ExposureCube::Aeff::value(cosTheta, phi);
}

} // namespace Likelihood
