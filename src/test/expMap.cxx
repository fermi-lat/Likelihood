/**
 * @file expMap.cxx
 * @brief Prototype standalone application for creating exposure maps used
 * by the Likelihood tool.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/test/expMap.cxx,v 1.3 2003/11/12 22:01:41 jchiang Exp $
 */

#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>

#include "optimizers/Exception.h"

#include "latResponse/IrfsFactory.h"

#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/ScData.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/RunParams.h"

using namespace Likelihood;

int main(int iargc, char* argv[]) {

#ifdef TRAP_FPE
   feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif

// Read in the command-line parameters using HOOPS
   std::string filename("expMap.par");
   delete argv[0];
   argv[0] = strdup(filename.c_str());

   RunParams params(iargc, argv);

// Set the region-of-interest.
   std::string roiCutsFile;
   params.getParam("ROI_cuts_file", roiCutsFile);
   RoiCuts::setCuts(roiCutsFile);

// Read in the pointing information.
   std::string scFile;
   params.getParam("Spacecraft_file", scFile);
   long scHdu;
   params.getParam("Spacecraft_file_hdu", scHdu);
   std::vector<std::string> scFiles;
   RunParams::resolve_fits_files(scFile, scFiles);
   std::vector<std::string>::const_iterator scIt = scFiles.begin();
   for ( ; scIt != scFiles.end(); scIt++) {
      ScData::readData(*scIt, scHdu);
   }

// Create the response functions.
   std::string responseFuncs;
   params.getParam("Response_functions", responseFuncs);
   latResponse::IrfsFactory irfsFactory;
   if (responseFuncs == "COMBINED_G25") {
      ResponseFunctions::addRespPtr(4, 
                                    irfsFactory.create("Glast25::Combined"));
   } else if (responseFuncs == "FRONT/BACK_G25") {
      ResponseFunctions::addRespPtr(2, irfsFactory.create("Glast25::Front"));
      ResponseFunctions::addRespPtr(3, irfsFactory.create("Glast25::Back"));
   } else if (responseFuncs == "TESTDC1") {
      ResponseFunctions::addRespPtr(1, irfsFactory.create("DC1::test"));
   } else if (responseFuncs == "FRONT") {
      ResponseFunctions::addRespPtr(5, irfsFactory.create("DC1::Front"));
   } else if (responseFuncs == "BACK") {
      ResponseFunctions::addRespPtr(6, irfsFactory.create("DC1::Back"));
   } else if (responseFuncs == "FRONT/BACK") {
      ResponseFunctions::addRespPtr(5, irfsFactory.create("DC1::Front"));
      ResponseFunctions::addRespPtr(6, irfsFactory.create("DC1::Back"));
   }

// Set the source region radius.  This should be larger than the
// radius of the region-of-interest.
   double sr_radius;
   params.getParam("Source_region_radius", sr_radius);
   RoiCuts *roiCuts = RoiCuts::instance();
   if (sr_radius < roiCuts->extractionRegion().radius() + 10.) {
      std::cerr << "The radius of the source region, " << sr_radius 
                << ", should be significantly larger (say by 10 deg) "
                << "than the ROI radius of " 
                << roiCuts->extractionRegion().radius() << std::endl;
      assert(sr_radius > roiCuts->extractionRegion().radius());
   }

// Get the other (hidden) map parameters.
   long nlong;
   params.getParam("number_of_longitude_points", nlong);
   long nlat;
   params.getParam("number_of_latitude_points", nlat);
   long nenergies;
   params.getParam("number_of_energies", nenergies);

// Create the exposure map.
   std::string exposureFile;
   params.getParam("Exposure_map_file", exposureFile);
   ExposureMap::computeMap(exposureFile, sr_radius, nlong, nlat, nenergies);

}
