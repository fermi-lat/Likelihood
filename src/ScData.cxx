/** 
 * @file ScData.cxx
 * @brief Implementation for the LAT spacecraft data class
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/ScData.cxx,v 1.13 2003/11/07 02:27:10 jchiang Exp $
 */

#include <vector>
#include <string>
#include <cmath>

#include "facilities/Util.h"

#include "Goodi/GoodiConstants.h"
#include "Goodi/DataIOServiceFactory.h"
#include "Goodi/DataFactory.h"
#include "Goodi/IDataIOService.h"
#include "Goodi/ISpacecraftData.h"
#include "Goodi/SCRow.h"

#include "astro/EarthCoordinate.h"

#include "latResponse/../src/Table.h"

#include "Likelihood/ScData.h"

namespace Likelihood {

// definitions of static data
std::vector<ScData::ScNtuple> ScData::vec;
std::string ScData::s_scFile = "";
int ScData::s_scHdu = 0;
ScData * ScData::s_instance = 0;
double ScData::s_tstep;

#ifdef USE_GOODI
void ScData::readData(std::string file, int hdu, bool clear) {
   facilities::Util::expandEnvVar(&file);

   s_scFile = file;
   s_scHdu = hdu;

   Goodi::DataIOServiceFactory iosvcCreator;
   Goodi::DataFactory dataCreator;
   Goodi::IDataIOService *ioService;
   ioService = iosvcCreator.create(file);
   Goodi::DataType datatype = Goodi::Spacecraft;
   Goodi::ISpacecraftData *scData = dynamic_cast<Goodi::ISpacecraftData *>
      (dataCreator.create(datatype, ioService));

   if (clear) vec.clear();
   bool done(false);
   while (!done) {
      const Goodi::SCRow *scRow = scData->nextSCrow(ioService, done);

      ScNtuple tuple;
      tuple.xAxis = astro::SkyDir(scRow->raSCX()*180./M_PI,
                                  scRow->decSCX()*180./M_PI);
      tuple.zAxis = astro::SkyDir(scRow->raSCZ()*180./M_PI,
                                  scRow->decSCZ()*180./M_PI);
      tuple.time = scRow->startTime();
      astro::EarthCoordinate earthCoord(scRow->latGeo()*180./M_PI,
                                        scRow->lonGeo()*180./M_PI);
      if (earthCoord.insideSAA()) {
         tuple.inSaa = 1;
      } else {
         tuple.inSaa = 0;
      }
      vec.push_back(tuple);
   }
   s_tstep = vec[1].time - vec[0].time;

   delete ioService;
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
#endif // USE_GOODI

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
