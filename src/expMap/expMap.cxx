/**
 * @file expMap.cxx
 * @brief Prototype standalone application for creating exposure maps used
 * by the Likelihood tool.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/expMap/expMap.cxx,v 1.7 2004/06/05 15:22:16 jchiang Exp $
 */

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <stdexcept>

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/Util.h"

using namespace Likelihood;

/**
 * @class ExpMap
 * @brief Class encapsulating methods for creating a Likelihood-specific
 * exposure map.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/expMap/expMap.cxx,v 1.7 2004/06/05 15:22:16 jchiang Exp $
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
   void setSourceRegion();
   void createExposureMap();
};

st_app::StAppFactory<ExpMap> myAppFactory;

ExpMap::ExpMap() : st_app::StApp(), m_helper(0), 
                   m_pars(st_app::StApp::getParGroup("expMap")) {
   try {
      m_pars.Prompt();
      m_pars.Save();
      m_helper = new AppHelpers(m_pars);
      ResponseFunctions::setEdispFlag(m_pars["use_energy_dispersion"]);
   } catch (std::exception &eObj) {
      std::cerr << eObj.what() << std::endl;
      std::exit(1);
   } catch (...) {
      std::cerr << "Caught unknown exception in ExpMap constructor." 
                << std::endl;
      std::exit(1);
   }
}

void ExpMap::run() {
   m_helper->setRoi();
   setSourceRegion();
   createExposureMap();
}

void ExpMap::setSourceRegion() {
   m_srRadius = m_pars["Source_region_radius"];
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
   long nlong = m_pars["number_of_longitude_points"];
   long nlat = m_pars["number_of_latitude_points"];
   long nenergies = m_pars["number_of_energies"];
// Exposure hypercube file.
   std::string expCubeFile = m_pars["exposure_cube_file"];
   if (expCubeFile != "none") {
      Util::file_ok(expCubeFile);
      PointSource::readExposureCube(expCubeFile);
   }
   std::string exposureFile = m_pars["Exposure_map_file"];
   ExposureMap::computeMap(exposureFile, m_srRadius, nlong, nlat, nenergies);
}
