/** 
 * @file ScData.cxx
 * @brief Implementation for the LAT spacecraft data class
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/ScData.cxx,v 1.19 2004/03/03 00:35:30 jchiang Exp $
 */

#include <cassert>
#include <vector>
#include <string>
#include <cmath>

#include "facilities/Util.h"

//#include "latResponse/../src/Table.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "astro/EarthCoordinate.h"

#include "Likelihood/ScData.h"

namespace Likelihood {

// definitions of static data
std::vector<ScData::ScNtuple> ScData::vec;
std::string ScData::s_scFile = "";
int ScData::s_scHdu = 0;
ScData * ScData::s_instance = 0;
double ScData::s_tstep;

#ifdef USE_TIP
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
#else
void ScData::readData(std::string file, int hdu, bool clear) {

   facilities::Util::expandEnvVar(&file);

   s_scFile = file;
   s_scHdu = hdu;

// read in the data (should check on file existence, etc., first...)
   latResponse::Table scTable;
   scTable.add_columns("SC_x0 SC_x1 SC_x2 SC_x SC_y SC_z time SAA_flag");
   scTable.read_FITS_table(file, hdu);

// repack into a more useful format
   if (clear) vec.clear();
   for (int i = 0; i < scTable[0].dim; i++) {
      ScNtuple tuple;

      tuple.xAxis = astro::SkyDir(Hep3Vector(scTable[0].val[i],
                                             scTable[1].val[i], 
                                             scTable[2].val[i]));
      tuple.zAxis = astro::SkyDir(Hep3Vector(scTable[3].val[i],
                                             scTable[4].val[i], 
                                             scTable[5].val[i]));
      tuple.time = scTable[6].val[i]; 
      tuple.inSaa = static_cast<int>(scTable[7].val[i]);

      vec.push_back(tuple);
   }

// Assume a constant time step.
   s_tstep = vec[1].time - vec[0].time;
}
#endif // USE_TIP

astro::SkyDir &ScData::zAxis(double time) {
   int indx = static_cast<int>((time - vec[0].time)/s_tstep);
   double frac = (time - vec[indx].time)/s_tstep;
   Hep3Vector zDir = frac*(vec[indx+1].zAxis.dir() - vec[indx].zAxis.dir())
      + vec[indx].zAxis.dir();
   m_zAxis = astro::SkyDir(zDir.unit());
   return m_zAxis;
}

astro::SkyDir &ScData::xAxis(double time) {
   int indx = static_cast<int>((time - vec[0].time)/s_tstep);
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

} // namespace Likelihood
