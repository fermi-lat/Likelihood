/**
 * @file expMap.cxx
 * @brief Prototype standalone application for creating exposure maps used
 * by the Likelihood tool.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/expMap/expMap.cxx,v 1.3 2004/04/06 01:17:00 jchiang Exp $
 */

#include <cmath>
#include <cstring>

#include <stdexcept>

#include "st_app/IApp.h"
#include "hoops/hoops_prompt_group.h"

#include "Likelihood/AppBase.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/PointSource.h"
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
 * $Header$
 */
class ExpMap : public AppBase {
public:
   ExpMap(hoops::IParGroup & pars) : AppBase(pars) {}
   virtual ~ExpMap() {}
   virtual void run();
private:
   double m_srRadius;
   void setSourceRegion();
   void createExposureMap();
};

/**
 * @class app
 * @brief Class (and object declaration) of boiler-plate code expected 
 * by st_app.
 */
class app : public st_app::IApp {
public:
   app() : st_app::IApp("expMap") {}
   virtual void run() {
      hoopsPrompt();
      hoopsSave();
      hoops::IParGroup & pars = hoopsGetParGroup();
      ExpMap expMapObject(pars);
      expMapObject.run();
   }
} my_app;

void ExpMap::run() {
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
