/** 
 * @file ScData.cxx
 * @brief Implementation for the LAT spacecraft data class
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/ScData.cxx,v 1.50 2009/05/15 04:10:41 jchiang Exp $
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

//#include "irfInterface/EfficiencyFactor.h"

#include "Likelihood/ScData.h"

namespace Likelihood {

void ScData::readData(std::string file, bool clear,
                      const std::string & sctable) {
   facilities::Util::expandEnvVar(&file);

   m_scFile = file;

   const tip::Table * scData = 
      tip::IFileSvc::instance().readTable(file, sctable);

   if (clear) {
      vec.clear();
   }

   double raSCX, decSCX;
   double raSCZ, decSCZ;
   tip::Table::ConstIterator it = scData->begin();
   tip::ConstTableRecord & scInterval = *it;
   for ( ; it != scData->end(); ++it) {
      ScNtuple tuple;
      scInterval["start"].get(tuple.time);
      scInterval["stop"].get(tuple.stoptime);
      scInterval["livetime"].get(tuple.livetime);
      scInterval["ra_scx"].get(raSCX);
      scInterval["dec_scx"].get(decSCX);
      tuple.xAxis = astro::SkyDir(raSCX, decSCX);
      scInterval["ra_scz"].get(raSCZ);
      scInterval["dec_scz"].get(decSCZ);
      tuple.zAxis = astro::SkyDir(raSCZ, decSCZ);
// Ensure that startTimes are contiguous.
      if (vec.size() > 1 && tuple.time < vec[vec.size()-2].time) {
         std::ostringstream message;
         message << "Likelihood::ScData: "
                 << "The start times in the spacecraft data are not "
                 << "contiguous.\n"
                 << "Previous time: " << vec[vec.size()-2].time << "\n"
                 << "Current time: " << tuple.time << "\n"
                 << "Current S/C file: " << m_scFile << "\n"
                 << "Check the ordering of your S/C files." 
                 << std::endl;
         throw std::runtime_error(message.str());
      }
      tuple.inSaa = 0;
      vec.push_back(tuple);
   }

   delete scData;
}

void ScData::readData(std::string file, double tstart, 
                      double tstop, bool clear,
                      const std::string & sctable) {
   static double maxIntervalSize(30);
   tstart -= 2*maxIntervalSize;
   tstop += 2*maxIntervalSize;
   facilities::Util::expandEnvVar(&file);

   m_scFile = file;

   std::ostringstream filter;
   filter << std::setprecision(10);
   filter << "(START >= " << tstart
          << ") && (STOP <= " << tstop << ")";
//    std::cout << "ScData::readData: using filter string: \n" 
//              << "  " << filter.str() << std::endl;
   const tip::Table * scData = 
      tip::IFileSvc::instance().readTable(file, sctable, filter.str());

   if (clear) {
      vec.clear();
//      irfInterface::EfficiencyFactor::clearFt2Data();
   }
//   irfInterface::EfficiencyFactor::readFt2File(file);

   double raSCX, decSCX;
   double raSCZ, decSCZ;
   tip::Table::ConstIterator it = scData->begin();
   tip::ConstTableRecord & scInterval = *it;
   for ( ; it != scData->end(); ++it) {
      ScNtuple tuple;
      scInterval["start"].get(tuple.time);
      scInterval["stop"].get(tuple.stoptime);
      scInterval["livetime"].get(tuple.livetime);
      scInterval["ra_scx"].get(raSCX);
      scInterval["dec_scx"].get(decSCX);
      tuple.xAxis = astro::SkyDir(raSCX, decSCX);
      scInterval["ra_scz"].get(raSCZ);
      scInterval["dec_scz"].get(decSCZ);
      tuple.zAxis = astro::SkyDir(raSCZ, decSCZ);
// Ensure that startTimes are contiguous.
      if (vec.size() > 1 && tuple.time < vec[vec.size()-2].time) {
         std::ostringstream message;
         message << "Likelihood::ScData: "
                 << "The start times in the spacecraft data are not "
                 << "contiguous.\n"
                 << "Previous time: " << vec[vec.size()-2].time << "\n"
                 << "Current time: " << tuple.time << "\n"
                 << "Current S/C file: " << m_scFile << "\n"
                 << "Check the ordering of your S/C files." 
                 << std::endl;
         throw std::runtime_error(message.str());
      }
      tuple.inSaa = 0;
      vec.push_back(tuple);
   }

   delete scData;
}

void ScData::readData(const std::vector<std::string> & scFiles, 
                      double tstart, double tstop,
                      const std::string & sctable) {
   vec.clear();
   std::vector<std::string>::const_iterator scIt = scFiles.begin();
   for ( ; scIt != scFiles.end(); scIt++) {
      st_facilities::Util::file_ok(*scIt);
      readData(*scIt, tstart, tstop, false, sctable);
   }
   if (vec.size() == 0) {
      throw std::runtime_error("No spacecraft time intervals were read in "
                               "for the desired range of FT1 data.");
   }
   for (size_t i(0); i < scFiles.size(); i++) {
   }
}

unsigned int ScData::time_index(double time) const {
   double tmin(vec.front().time);
   double tmax(vec.back().stoptime);
   double tol(1e-5);
   if (time < tmin - tol || time > tmax + tol) {
      std::ostringstream message;
      message << "Requested time of " << time << " "
              << "lies outside the range of valid times in the "
              << "pointing/livetime history: " 
              << tmin << " to " << tmax << "MET s";
      throw std::runtime_error(message.str());
   }
   ScNtuple my_vec;
   my_vec.time = time;
   std::vector<ScNtuple>::const_iterator it
      = std::upper_bound(vec.begin(), vec.end(), my_vec, less_than_time);
   unsigned int indx = it - vec.begin() - 1;
   return indx;
}

astro::SkyDir ScData::zAxis(double time) const {
   size_t indx = time_index(time);
   indx = std::min(indx, vec.size()-2);
   double frac = (time - vec[indx].time)
      /(vec.at(indx+1).time - vec.at(indx).time);
   Hep3Vector zDir = frac*(vec[indx+1].zAxis.dir() - vec[indx].zAxis.dir())
      + vec[indx].zAxis.dir();
   return astro::SkyDir(zDir.unit());
}

astro::SkyDir ScData::xAxis(double time) const {
   size_t indx = time_index(time);
   indx = std::min(indx, vec.size()-2);
   double frac = (time - vec[indx].time)
      /(vec.at(indx+1).time - vec.at(indx).time);
   Hep3Vector xDir = frac*(vec[indx+1].xAxis.dir() - vec[indx].xAxis.dir())
      + vec[indx].xAxis.dir();
   return astro::SkyDir(xDir.unit());
}

std::pair<ScData::Iterator, ScData::Iterator> 
ScData::bracketInterval(double startTime, double stopTime) const {
   ScNtuple startTuple;
   startTuple.time = startTime;
   ScData::Iterator lowerBound 
      = std::lower_bound(vec.begin(), vec.end(), startTuple,
                         Likelihood::ScData::less_than_time);
   if (lowerBound == vec.end()) {
      std::ostringstream message;
      message << "Likelihood::ScData::bracketInterval:\nStart time " 
              << startTime << " is out-of-range for "
              << "existing spacecraft data time range: (" 
              << (*vec.begin()).time
              << ", " << (*(vec.end()-1)).time << ")";
      throw std::out_of_range(message.str());
   }
   if ((*lowerBound).time != startTime) --lowerBound;

   ScNtuple stopTuple;
   stopTuple.time = stopTime;
   ScData::Iterator upperBound 
      = std::upper_bound(vec.begin(), vec.end(), stopTuple,
                         Likelihood::ScData::less_than_time);
   if (upperBound == vec.end()) {
      std::ostringstream message;
      message << "Likelihood::ScData::bracketInterval:\nStop time " 
              << stopTime << " is out-of-range for "
              << "existing spacecraft data time range: (" 
              << (*vec.begin()).time
              << ", " << (*(vec.end()-1)).time << ")";
      throw std::out_of_range(message.str());
   }
   return std::make_pair(lowerBound, upperBound);
}

bool ScData::less_than_time(const ScNtuple & scDatum1, 
                            const ScNtuple & scDatum2) {
   return scDatum1.time < scDatum2.time;
}

} // namespace Likelihood
