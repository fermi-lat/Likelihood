/** 
 * @file ScData.cxx
 * @brief Implementation for the LAT spacecraft data class
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/ScData.cxx,v 1.9 2003/10/22 04:30:33 jchiang Exp $
 */

#include <vector>
#include <string>
#include <cmath>

#include "latResponse/../src/Table.h"

#include "Likelihood/ScData.h"

namespace Likelihood {

// definitions of static data
std::vector<ScData::ScNtuple> ScData::vec;
std::string ScData::s_scFile = "";
int ScData::s_scHdu = 0;
ScData * ScData::s_instance = 0;
double ScData::s_tstep;

void ScData::readData(const std::string &file, int hdu) {
   s_scFile = file;
   s_scHdu = hdu;

// read in the data (should check on file existence, etc., first...)
   latResponse::Table scTable;
   scTable.add_columns("SC_x0 SC_x1 SC_x2 SC_x SC_y SC_z time SAA_flag");
   scTable.read_FITS_table(file, hdu);

// repack into a more useful format
//    vec.clear();
//    vec.reserve(scTable[0].dim);
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

astro::SkyDir &ScData::zAxis(double time) {
   int indx = static_cast<int>((time - vec[0].time)/s_tstep);
   double frac = (time - vec[indx].time)/s_tstep;
   Hep3Vector zDir = frac*(vec[indx+1].zAxis.dir() - vec[indx].zAxis.dir())
      + vec[indx].zAxis.dir();
   m_zAxis = astro::SkyDir(zDir);
   return m_zAxis;
}

astro::SkyDir &ScData::xAxis(double time) {
   int indx = static_cast<int>((time - vec[0].time)/s_tstep);
   double frac = (time - vec[indx].time)/s_tstep;
   Hep3Vector xDir = frac*(vec[indx+1].xAxis.dir() - vec[indx].xAxis.dir())
      + vec[indx].xAxis.dir();
   m_xAxis = astro::SkyDir(xDir);
   return m_xAxis;
}

ScData * ScData::instance() {
   if (s_instance == 0) {
      s_instance = new ScData();
   }
   return s_instance;
}

} // namespace Likelihood
