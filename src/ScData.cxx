/** 
 * @file ScData.cxx
 * @brief Implementation for the LAT spacecraft data class
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/ScData.cxx,v 1.54 2009/06/02 20:51:39 jchiang Exp $
 */

#include <cmath>

#include <algorithm>
#include <iomanip>
#include <string>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "facilities/Util.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "astro/EarthCoordinate.h"

#include "st_facilities/Util.h"

#include "irfInterface/EfficiencyFactor.h"

#include "Likelihood/ScData.h"

namespace Likelihood {

void ScData::readData(std::string scfile, double tstart, double tstop,
                      bool clear, const std::string & sctable) {
   static double maxIntervalSize(30);
   tstart -= 2*maxIntervalSize;
   tstop += 2*maxIntervalSize;
   facilities::Util::expandEnvVar(&scfile);

   std::ostringstream filter;
   filter << std::setprecision(10);
   filter << "(START >= " << tstart
          << ") && (STOP <= " << tstop << ")";

   const tip::Table * scData = 
      tip::IFileSvc::instance().readTable(scfile, sctable, filter.str());

   if (clear) {
      clear_arrays();
   }

   double start, stop, livetime;
   double raSCX, decSCX;
   double raSCZ, decSCZ;
   tip::Table::ConstIterator it = scData->begin();
   tip::ConstTableRecord & scInterval = *it;
   for ( ; it != scData->end(); ++it) {
      scInterval["start"].get(start);
      scInterval["stop"].get(stop);
      scInterval["livetime"].get(livetime);
      scInterval["ra_scx"].get(raSCX);
      scInterval["dec_scx"].get(decSCX);
      scInterval["ra_scz"].get(raSCZ);
      scInterval["dec_scz"].get(decSCZ);
// Ensure that start times are monotonically increasing.
      if (m_start.size() > 1 && start < m_start.back()) {
         std::ostringstream message;
         message << "Likelihood::ScData: "
                 << "The start times in the spacecraft data are not "
                 << "monitonically increasing.\n"
                 << "Previous time: " << m_start.back() << "\n"
                 << "Current time: " << start << "\n"
                 << "Current S/C file: " << scfile << "\n"
                 << "Check the ordering of your S/C files." 
                 << std::endl;
         throw std::runtime_error(message.str());
      }
      m_start.push_back(start);
      m_stop.push_back(stop);
      m_livetime.push_back(livetime);
      m_xAxis.push_back(astro::SkyDir(raSCX, decSCX));
      m_zAxis.push_back(astro::SkyDir(raSCZ, decSCZ));
   }
   delete scData;
}

void ScData::readData(const std::vector<std::string> & scFiles, 
                      double tstart, double tstop,
                      const std::string & sctable) {
   clear_arrays();
   std::vector<std::string>::const_iterator scIt = scFiles.begin();
   for ( ; scIt != scFiles.end(); scIt++) {
      st_facilities::Util::file_ok(*scIt);
      readData(*scIt, tstart, tstop, false, sctable);
   }
   if (m_start.size() == 0) {
      throw std::runtime_error("No spacecraft time intervals were read in "
                               "for the desired range of FT1 data.");
   }
}

size_t ScData::time_index(double time) const {
   double tmin(m_start.front());
   double tmax(m_stop.back());
   double tol(1e-5);
   if (time < tmin - tol || time > tmax + tol) {
      std::ostringstream message;
      message << "Requested time of " << time << " "
              << "lies outside the range of valid times in the "
              << "pointing/livetime history: " 
              << tmin << " to " << tmax << "MET s";
      throw std::runtime_error(message.str());
   }
   std::vector<double>::const_iterator it 
      = std::upper_bound(m_start.begin(), m_start.end(), time);
   size_t indx = it - m_start.begin() - 1;
   return indx;
}

double ScData::livetimefrac(double time) const {
   size_t indx = time_index(time);
   return m_livetime.at(indx)/(m_stop.at(indx) - m_start.at(indx));
}

astro::SkyDir ScData::xAxis(double time) const {
   size_t indx = time_index(time);
   indx = std::min(indx, m_xAxis.size() - 2);
   double frac = (time - start(indx))/(start(indx+1) - start(indx));
   Hep3Vector zDir = frac*(xAxis(indx+1).dir() - xAxis(indx).dir())
      + xAxis(indx).dir();
   return astro::SkyDir(zDir.unit());
}

astro::SkyDir ScData::zAxis(double time) const {
   size_t indx = time_index(time);
   indx = std::min(indx, m_zAxis.size() - 2);
   double frac = (time - start(indx))/(start(indx+1) - start(indx));
   Hep3Vector zDir = frac*(zAxis(indx+1).dir() - zAxis(indx).dir())
      + zAxis(indx).dir();
   return astro::SkyDir(zDir.unit());
}

void ScData::clear_arrays() {
   m_start.clear();
   m_stop.clear();
   m_livetime.clear();
   m_xAxis.clear();
   m_zAxis.clear();
}

} // namespace Likelihood
