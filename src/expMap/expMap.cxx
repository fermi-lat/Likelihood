/**
 * @file expMap.cxx
 * @brief Prototype standalone application for creating exposure maps used
 * by the Likelihood tool.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/expMap/expMap.cxx,v 1.18 2005/01/19 02:15:01 jchiang Exp $
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
#include "Likelihood/ExposureMap.h"
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/expMap/expMap.cxx,v 1.18 2005/01/19 02:15:01 jchiang Exp $
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
private:
   AppHelpers * m_helper;
   st_app::AppParGroup & m_pars;
   double m_srRadius;
   void promptForParameters();
   void setSourceRegion();
   void createExposureMap();
};

st_app::StAppFactory<ExpMap> myAppFactory;

ExpMap::ExpMap() : st_app::StApp(), m_helper(0), 
                   m_pars(st_app::StApp::getParGroup("expMap")) {}

void ExpMap::run() {
   promptForParameters();
   Likelihood::Verbosity::instance(m_pars["chatter"]);
   m_helper = new AppHelpers(m_pars);
   m_helper->readScData();
   ResponseFunctions::setEdispFlag(m_pars["use_energy_dispersion"]);
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
      AppHelpers::checkTimeCuts(m_pars["evfile"], "EVENTS",
                                m_pars["exposure_cube_file"], "");
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
   RoiCuts *roiCuts = RoiCuts::instance();
   if (Likelihood::print_output() &&
       m_srRadius < roiCuts->extractionRegion().radius() + 10.) {
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
   long nlong = m_pars["number_of_longitude_points"];
   long nlat = m_pars["number_of_latitude_points"];
   long nenergies = m_pars["number_of_energies"];
// Exposure hypercube file.
   std::string expCubeFile = m_pars["exposure_cube_file"];
   if (expCubeFile != "none") {
      st_facilities::Util::file_ok(expCubeFile);
      ExposureCube::readExposureCube(expCubeFile);
   }
   std::string exposureFile = m_pars["outfile"];
   ExposureMap::computeMap(exposureFile, m_srRadius, nlong, nlat, nenergies);

   std::auto_ptr<tip::Image> 
      image(tip::IFileSvc::instance().editImage(exposureFile, ""));
   Likelihood::RoiCuts::instance()->writeDssKeywords(image->getHeader());
   Likelihood::RoiCuts::instance()->writeGtiExtension(exposureFile);
}
