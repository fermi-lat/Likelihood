/**
 * @file expMap.cxx
 * @brief Prototype standalone application for creating exposure maps used
 * by the Likelihood tool.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/test/expMap.cxx,v 1.6 2004/01/10 22:53:10 jchiang Exp $
 */

#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>

#include "facilities/Util.h"

#include "hoops/hoops_exception.h"
#include "hoopsUtil/RunParams.h"

#include "optimizers/Exception.h"

#include "latResponse/IrfsFactory.h"

#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/ScData.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/Exception.h"

using namespace Likelihood;
using latResponse::irfsFactory;

namespace {
   bool fileExists(const std::string &filename) {
      std::ifstream file(filename.c_str());
      return file.is_open();
   }

   void file_ok(std::string filename) {
      facilities::Util::expandEnvVar(&filename);
      if (::fileExists(filename)) {
         return;
      } else {
         std::cout << "likelihood::main:\n"
                   << "File not found: " << filename
                   << std::endl;
         exit(-1);
      }
   }

   void readLines(std::string inputFile, 
                  std::vector<std::string> &lines) {
      
      facilities::Util::expandEnvVar(&inputFile);
      
      std::ifstream file(inputFile.c_str());
      lines.clear();
      std::string line;
      while (std::getline(file, line, '\n')) {
         if (line != "" && line != " ") { //skip (most) blank lines
            lines.push_back(line);
         }
      }
   }

   void resolve_fits_files(std::string filename, 
                           std::vector<std::string> &files) {
      
      facilities::Util::expandEnvVar(&filename);
      files.clear();
      
// Read the first line of the file and see if the first 6 characters
// are "SIMPLE".  If so, then we assume it's a FITS file.
      std::ifstream file(filename.c_str());
      std::string firstLine;
      std::getline(file, firstLine, '\n');
      if (firstLine.find("SIMPLE") == 0) {
// This is a FITS file. Return that as the sole element in the files
// vector.
         files.push_back(filename);
         return;
      } else {
// filename contains a list of fits files.
         ::readLines(filename, files);
         return;
      }
   }

}

int main(int iargc, char* argv[]) {

#ifdef TRAP_FPE
   feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif

// Read in the command-line parameters using HOOPS
   strcpy(argv[0], "expMap");

   try {
      hoopsUtil::RunParams params(iargc, argv);

// Set the region-of-interest.
      std::string roiCutsFile;
      params.getParam("ROI_cuts_file", roiCutsFile);
      ::file_ok(roiCutsFile);
      RoiCuts::setCuts(roiCutsFile);

// Read in the pointing information.
      std::string scFile;
      params.getParam("Spacecraft_file", scFile);
      long scHdu;
      params.getParam("Spacecraft_file_hdu", scHdu);
      std::vector<std::string> scFiles;
      ::file_ok(scFile);
      ::resolve_fits_files(scFile, scFiles);
      std::vector<std::string>::const_iterator scIt = scFiles.begin();
      for ( ; scIt != scFiles.end(); scIt++) {
         ::file_ok(*scIt);
         ScData::readData(*scIt, scHdu);
      }

// Create the response functions.
      std::string responseFuncs;
      params.getParam("Response_functions", responseFuncs);

      std::map< std::string, std::vector<std::string> > responseIds;
      responseIds["FRONT"].push_back("DC1::Front");
      responseIds["BACK"].push_back("DC1::Back");
      responseIds["FRONT/BACK"].push_back("DC1::Front");
      responseIds["FRONT/BACK"].push_back("DC1::Back");

      if (responseIds.count(responseFuncs)) {
         std::vector<std::string> &resps = responseIds[responseFuncs];
         for (unsigned int i = 0; i < resps.size(); i++) {
            ResponseFunctions::addRespPtr(i, irfsFactory().create(resps[i]));
         }
      } else {
         std::cerr << "Invalid response function choice: "
                   << responseFuncs << std::endl;
         exit(-1);
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

   } catch (hoops::Hexception &eObj) {
      std::cout << "HOOPS exception code: "
                << eObj.Msg() << std::endl;
   } catch (Exception &eObj) {
      std::cout << eObj.what() << std::endl;
   }
}
