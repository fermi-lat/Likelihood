#include <vector>
#include <string>
#include <cmath>

#include "../Likelihood/ScData.h"
#include "../Likelihood/Table.h"

namespace Likelihood {

// definitions of static data
std::vector<ScData::ScNtuple> ScData::vec;
std::string ScData::m_scFile = "";
int ScData::m_scHdu = 0;
ScData * ScData::s_instance = 0;

void ScData::readData(const std::string &file, int hdu) {
   m_scFile = file;
   m_scHdu = hdu;

// read in the data (should check on file existence, etc., first...)
   Table scTable;
   scTable.add_columns("SC_x0 SC_x1 SC_x2 SC_x SC_y SC_z time SAA_flag");
   scTable.read_FITS_table(file, hdu);

// repack into a more useful format
   vec.reserve(scTable[0].dim);
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
   }
}

ScData * ScData::instance() {
   if (s_instance == 0) {
      s_instance = new ScData();
   }
   return s_instance;
}

} // namespace Likelihood
