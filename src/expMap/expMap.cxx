/**
 * @file expMap.cxx
 * @brief Prototype standalone application for creating exposure maps used
 * by the Likelihood tool.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/expMap/expMap.cxx,v 1.30 2005/04/21 19:00:26 jchiang Exp $
 */

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <memory>
#include <stdexcept>

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "st_facilities/Util.h"

#include "tip/IFileSvc.h"
#include "tip/Image.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/ExposureCube.h"
#include "Likelihood/Observation.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/RoiCuts.h"

#include "Verbosity.h"

using namespace Likelihood;

/**
 * @class ExpMap
 * @brief Class encapsulating methods for creating a Likelihood-specific
 * exposure map.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/expMap/expMap.cxx,v 1.30 2005/04/21 19:00:26 jchiang Exp $
 */
class ExpMap : public st_app::StApp {
public:
   ExpMap();
   virtual ~ExpMap() throw() {
      try {
         delete m_helper;
      } catch (std::exception &eObj) {
         std::cerr << eObj.what() << std::endl;
      } catch (...) {
      }
   }
   virtual void run();
   virtual void banner() const {}
private:
   AppHelpers * m_helper;
   st_app::AppParGroup & m_pars;
   double m_srRadius;
   void promptForParameters();
   void setSourceRegion();
   void createExposureMap();
};

st_app::StAppFactory<ExpMap> myAppFactory("gtexpmap");

ExpMap::ExpMap() : st_app::StApp(), m_helper(0), 
                   m_pars(st_app::StApp::getParGroup("gtexpmap")) {}

void ExpMap::run() {
   promptForParameters();
   Likelihood::Verbosity::instance(m_pars["chatter"]);
   m_helper = new AppHelpers(&m_pars);
   m_helper->readScData();
   bool useEdisp = m_pars["use_energy_dispersion"];
   ResponseFunctions & respFuncs =
      const_cast<ResponseFunctions &>(m_helper->observation().respFuncs());
   respFuncs.setEdispFlag(useEdisp);
   m_helper->setRoi();
   setSourceRegion();
   createExposureMap();
}

void ExpMap::promptForParameters() {
   m_pars.Prompt("evfile");
   m_pars.Prompt("scfile");
   m_pars.Prompt("exposure_cube_file");
   std::string expCubeFile = m_pars["exposure_cube_file"];
   if (expCubeFile != "none") {
      std::vector<std::string> eventFiles;
      st_facilities::Util::resolve_fits_files(m_pars["evfile"], eventFiles);
      AppHelpers::checkTimeCuts(eventFiles, "EVENTS",
                                m_pars["exposure_cube_file"], "Exposure");
   }
   m_pars.Prompt("outfile");
   AppHelpers::checkOutputFile(m_pars["clobber"], m_pars["outfile"]);
   m_pars.Prompt("rspfunc");
   m_pars.Prompt("source_region_radius");
   m_pars.Prompt("number_of_longitude_points");
   m_pars.Prompt("number_of_latitude_points");
   m_pars.Prompt("number_of_energies");
   m_pars.Save();
}

void ExpMap::setSourceRegion() {
   m_srRadius = m_pars["source_region_radius"];
   const RoiCuts & roiCuts = m_helper->observation().roiCuts();
   if (Likelihood::print_output() &&
       m_srRadius < roiCuts.extractionRegion().radius() + 10.) {
      std::cerr << "The radius of the source region, " << m_srRadius 
                << ", should be significantly larger (say by 10 deg) "
                << "than the ROI radius of " 
                << roiCuts.extractionRegion().radius() << std::endl;
      if (m_srRadius < roiCuts.extractionRegion().radius()) {
         std::ostringstream message;
         message << "The source region radius, " << m_srRadius 
                 << ", should be larger than the ROI radius, "
                 << roiCuts.extractionRegion().radius();
         throw std::out_of_range(message.str());
      }
   }
}

void ExpMap::createExposureMap() {
   long nlong = m_pars["number_of_longitude_points"];
   long nlat = m_pars["number_of_latitude_points"];
   long nenergies = m_pars["number_of_energies"];
// Exposure hypercube file.
   std::string expCubeFile = m_pars["exposure_cube_file"];
   if (expCubeFile != "none") {
      st_facilities::Util::file_ok(expCubeFile);
      ExposureCube & expCube = 
         const_cast<ExposureCube &>(m_helper->observation().expCube());
      expCube.readExposureCube(expCubeFile);
   }
   std::string exposureFile = m_pars["outfile"];
   const Observation & observation = m_helper->observation();
   const RoiCuts & roiCuts = observation.roiCuts();
   m_helper->observation().expMap().computeMap(exposureFile, observation,
                                               m_srRadius, nlong, nlat,
                                               nenergies); 

   std::auto_ptr<tip::Image> 
      image(tip::IFileSvc::instance().editImage(exposureFile, ""));
   roiCuts.writeDssKeywords(image->getHeader());
   roiCuts.writeGtiExtension(exposureFile);
}
