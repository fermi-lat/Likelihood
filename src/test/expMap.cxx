/**
 * @file expMap.cxx
 * @brief Prototype standalone application for creating exposure maps used
 * by the Likelihood tool.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/test/likelihood.cxx,v 1.5 2003/11/07 02:27:11 jchiang Exp $
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
   std::string roiCutsFile = params.string_par("ROI_cuts_file");
   RoiCuts::setCuts(roiCutsFile);

// Read in the pointing information.
   std::string scFile = params.string_par("Spacecraft_file");
   int scHdu = static_cast<int>(params.long_par("Spacecraft_file_hdu"));
   ScData::readData(scFile, scHdu);

// Create the response functions.
   std::string responseFuncs = params.string_par("Response_functions");
   latResponse::IrfsFactory irfsFactory;
   if (responseFuncs == "COMBINED") {
      ResponseFunctions::addRespPtr(4, 
                                    irfsFactory.create("Glast25::Combined"));
   } else if (responseFuncs == "FRONT/BACK") {
      ResponseFunctions::addRespPtr(2, irfsFactory.create("Glast25::Front"));
      ResponseFunctions::addRespPtr(3, irfsFactory.create("Glast25::Back"));
   }

// Set the source region radius.  This should be larger than the
// radius of the region-of-interest.
   double sr_radius = params.double_par("Source_region_radius");
   RoiCuts *roiCuts = RoiCuts::instance();
   if (sr_radius < roiCuts->extractionRegion().radius() + 10.) {
      std::cerr << "The radius of the source region, " << sr_radius 
                << ", should be significantly larger (say 10 deg) "
                << "than the ROI radius of " 
                << roiCuts->extractionRegion().radius() << std::endl;
      assert(sr_radius > roiCuts->extractionRegion().radius());
   }

// Get the other (hidden) map parameters.
   int nlong = params.long_par("number_of_longitude_points");
   int nlat = params.long_par("number_of_latitude_points");
   int nenergies = params.long_par("number_of_energies");

// Create the exposure map.
   std::string exposureFile = params.string_par("Exposure_map_file");
   ExposureMap::computeMap(exposureFile, sr_radius, nlong, nlat, nenergies);

}
