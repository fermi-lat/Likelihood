/**
 * @file expMap.cxx
 * @brief Prototype standalone application for creating exposure maps used
 * by the Likelihood tool.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/expMap/expMap.cxx,v 1.2 2004/04/05 22:01:56 jchiang Exp $
 */

#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <cmath>
#include <cstring>

#include <stdexcept>

#include "facilities/Util.h"
#include "hoops/hoops_prompt_group.h"
#include "st_app/IApp.h"

#include "latResponse/IrfsFactory.h"

#include "Likelihood/PointSource.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/ScData.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/Util.h"

using namespace Likelihood;
using latResponse::irfsFactory;

class ExpMap : public st_app::IApp {
public:
   ExpMap() : st_app::IApp("expMap") {}
   virtual ~ExpMap() throw() {}
   virtual void run();

private:
   double m_srRadius;
   void setUp();
   void tearDown();
   void setRoi();
   void readScData();
   void createResponseFunctions();
   void setSourceRegion();
   void createExposureMap();
};

st_app::IApp * my_application = new ExpMap();

void ExpMap::run() {
   setUp();
   setRoi();
   readScData();
   createResponseFunctions();
   setSourceRegion();
   createExposureMap();
   tearDown();
}

void ExpMap::setUp() {
#ifdef TRAP_FPE
   feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif
   hoopsPrompt();
   hoopsSave();
}

void ExpMap::tearDown() {
}

void ExpMap::setRoi() {
   hoops::IParGroup & pars = hoopsGetParGroup();
   std::string roiCutsFile = pars["ROI_cuts_file"];
   Util::file_ok(roiCutsFile);
   RoiCuts::setCuts(roiCutsFile);
}

void ExpMap::readScData() {
   hoops::IParGroup & pars = hoopsGetParGroup();
   std::string scFile = pars["Spacecraft_file"];
   long scHdu = pars["Spacecraft_file_hdu"];
   std::vector<std::string> scFiles;
   Util::file_ok(scFile);
   Util::resolve_fits_files(scFile, scFiles);
   std::vector<std::string>::const_iterator scIt = scFiles.begin();
   for ( ; scIt != scFiles.end(); scIt++) {
      Util::file_ok(*scIt);
      ScData::readData(*scIt, scHdu);
   }
}

void ExpMap::createResponseFunctions() {
   hoops::IParGroup & pars = hoopsGetParGroup();
   std::string responseFuncs = pars["Response_functions"];
   std::map< std::string, std::vector<std::string> > responseIds;
   responseIds["FRONT"].push_back("DC1::Front");
   responseIds["BACK"].push_back("DC1::Back");
   responseIds["FRONT/BACK"].push_back("DC1::Front");
   responseIds["FRONT/BACK"].push_back("DC1::Back");
   responseIds["GLAST25"].push_back("Glast25::Front");
   responseIds["GLAST25"].push_back("Glast25::Back");
   if (responseIds.count(responseFuncs)) {
      std::vector<std::string> &resps = responseIds[responseFuncs];
      for (unsigned int i = 0; i < resps.size(); i++) {
         ResponseFunctions::addRespPtr(i, irfsFactory().create(resps[i]));
      }
   } else {
      throw std::invalid_argument("Invalid response function choice: "
                                  + responseFuncs);
   }
}

void ExpMap::setSourceRegion() {
   hoops::IParGroup & pars = hoopsGetParGroup();
   m_srRadius = pars["Source_region_radius"];
   RoiCuts *roiCuts = RoiCuts::instance();
   if (m_srRadius < roiCuts->extractionRegion().radius() + 10.) {
      std::cerr << "The radius of the source region, " << m_srRadius 
                << ", should be significantly larger (say by 10 deg) "
                << "than the ROI radius of " 
                << roiCuts->extractionRegion().radius() << std::endl;
      if (m_srRadius < roiCuts->extractionRegion().radius()) {
         std::ostringstream message;
         message << "The source region radius, " << m_srRadius 
                 << ", should be larger than the ROI radius, "
                 << roiCuts->extractionRegion().radius();
         throw std::out_of_range(message.str());
      }
   }
}

void ExpMap::createExposureMap() {
   hoops::IParGroup & pars = hoopsGetParGroup();
   long nlong = pars["number_of_longitude_points"];
   long nlat = pars["number_of_latitude_points"];
   long nenergies = pars["number_of_energies"];
// Exposure hypercube file.
   std::string expCubeFile = pars["exposure_cube_file"];
   if (expCubeFile != "none") {
      Util::file_ok(expCubeFile);
      PointSource::readExposureCube(expCubeFile);
   }
   std::string exposureFile = pars["Exposure_map_file"];
   ExposureMap::computeMap(exposureFile, m_srRadius, nlong, nlat, nenergies);
}
