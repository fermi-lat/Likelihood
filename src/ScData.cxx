/** 
 * @file ScData.cxx
 * @brief Implementation for the LAT spacecraft data class
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/ScData.cxx,v 1.26 2004/07/21 04:00:13 jchiang Exp $
 */

#include <cassert>
#include <cmath>

#include <algorithm>
#include <string>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "facilities/Util.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "irfUtil/Util.h"

#include "astro/EarthCoordinate.h"

#include "Likelihood/ScData.h"

namespace Likelihood {

// definitions of static data
std::vector<ScData::ScNtuple> ScData::vec;
std::string ScData::s_scFile = "";
int ScData::s_scHdu = 0;
ScData * ScData::s_instance = 0;
double ScData::s_tstep;

void ScData::readData(std::string file, int hdu, bool clear) {
   facilities::Util::expandEnvVar(&file);

   s_scFile = file;
   s_scHdu = hdu;

   tip::Table * scData = tip::IFileSvc::instance().editTable(file, "ext1");

   if (clear) vec.clear();

   double raSCX, decSCX;
   double raSCZ, decSCZ;
//   double lonGeo, latGeo;
   tip::Table::Iterator it = scData->begin();
   tip::Table::Record & scInterval = *it;
   for ( ; it != scData->end(); ++it) {
      ScNtuple tuple;
      scInterval["start"].get(tuple.time);
      scInterval["ra_scx"].get(raSCX);
      scInterval["dec_scx"].get(decSCX);
      tuple.xAxis = astro::SkyDir(raSCX, decSCX);
      scInterval["ra_scz"].get(raSCZ);
      scInterval["dec_scz"].get(decSCZ);
      tuple.zAxis = astro::SkyDir(raSCZ, decSCZ);
// Ensure that startTimes are contiguous.
      if (vec.size() > 1 && tuple.time < vec[vec.size()-2].time) {
         std::cerr << "Likelihood::ScData: "
                   << "The start times in the spacecraft data are not "
                   << "contiguous.\n"
                   << "Previous time: " << vec[vec.size()-2].time << "\n"
                   << "Current time: " << tuple.time << "\n"
                   << "Current S/C file: " << s_scFile << "\n"
                   << "Check the ordering of your S/C files." 
                   << std::endl;
         assert(tuple.time > vec[vec.size()-2].time);
      }
//       scInterval["lon_geo"].get(lonGeo);
//       scInterval["lat_geo"].get(latGeo);
//       astro::EarthCoordinate earthCoord(latGeo, lonGeo);
//       if (earthCoord.insideSAA()) {
//          tuple.inSaa = 1;
//       } else {
//          tuple.inSaa = 0;
//       }
      tuple.inSaa = 0;
      vec.push_back(tuple);
   }
   s_tstep = vec[1].time - vec[0].time;

   delete scData;
}         

astro::SkyDir &ScData::zAxis(double time) {
   int indx = static_cast<int>((time - vec[0].time)/s_tstep);
   indx = std::min(static_cast<unsigned int>(indx), vec.size()-2);
   double frac = (time - vec[indx].time)/s_tstep;
   Hep3Vector zDir = frac*(vec[indx+1].zAxis.dir() - vec[indx].zAxis.dir())
      + vec[indx].zAxis.dir();
   m_zAxis = astro::SkyDir(zDir.unit());
   return m_zAxis;
}

astro::SkyDir &ScData::xAxis(double time) {
   int indx = static_cast<int>((time - vec[0].time)/s_tstep);
   indx = std::min(static_cast<unsigned int>(indx), vec.size()-2);
   double frac = (time - vec[indx].time)/s_tstep;
   Hep3Vector xDir = frac*(vec[indx+1].xAxis.dir() - vec[indx].xAxis.dir())
      + vec[indx].xAxis.dir();
   m_xAxis = astro::SkyDir(xDir.unit());
   return m_xAxis;
}

ScData * ScData::instance() {
   if (s_instance == 0) {
      s_instance = new ScData();
   }
   return s_instance;
}

std::pair<ScData::Iterator, ScData::Iterator> 
ScData::bracketInterval(double startTime, double stopTime) {
   
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
