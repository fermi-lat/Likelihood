/**
 * @file expMap.cxx
 * @brief Prototype standalone application for creating exposure maps used
 * by the Likelihood tool.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/expMap/expMap.cxx,v 1.42 2008/03/24 22:45:44 jchiang Exp $
 */

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <memory>
#include <stdexcept>

#include "st_stream/StreamFormatter.h"

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

using namespace Likelihood;

/**
 * @class ExpMap
 * @brief Class encapsulating methods for creating a Likelihood-specific
 * exposure map.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/expMap/expMap.cxx,v 1.42 2008/03/24 22:45:44 jchiang Exp $
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
   virtual void banner() const;
private:
   AppHelpers * m_helper;
   st_app::AppParGroup & m_pars;
   double m_srRadius;
   void promptForParameters();
   void setSourceRegion();
   void createExposureMap();
   static std::string s_cvs_id;
};

st_app::StAppFactory<ExpMap> myAppFactory("gtexpmap");

std::string ExpMap::s_cvs_id("$Name:  $");

ExpMap::ExpMap() : st_app::StApp(), m_helper(0), 
                   m_pars(st_app::StApp::getParGroup("gtexpmap")) {
   setVersion(s_cvs_id);
   m_pars.setSwitch("submap");
   m_pars.setCase("submap", "yes", "nlongmin");
   m_pars.setCase("submap", "yes", "nlongmax");
   m_pars.setCase("submap", "yes", "nlatmin");
   m_pars.setCase("submap", "yes", "nlatmax");
}

void ExpMap::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

void ExpMap::run() {
   st_stream::StreamFormatter formatter("gtexpmap", "run", 2);
   formatter.warn() << "The exposure maps generated by this tool are meant\n"
                    << "to be used for *unbinned* likelihood analysis only.\n"
                    << "Do not use them for binned analyses." << std::endl;
   promptForParameters();
   m_helper = new AppHelpers(&m_pars, "UNBINNED");
//   bool useEdisp = m_pars["edisp"];
   bool useEdisp = false;
   ResponseFunctions & respFuncs =
      const_cast<ResponseFunctions &>(m_helper->observation().respFuncs());
   respFuncs.setEdispFlag(useEdisp);
   m_helper->setRoi();
   m_helper->readScData();
   setSourceRegion();
   createExposureMap();
}

void ExpMap::promptForParameters() {
   m_pars.Prompt("evfile");
   m_pars.Prompt("scfile");
   m_pars.Prompt("expcube");
   std::string expCubeFile = m_pars["expcube"];
   std::string evtable = m_pars["evtable"];
   if (expCubeFile != "none") {
      std::vector<std::string> eventFiles;
      st_facilities::Util::resolve_fits_files(m_pars["evfile"], eventFiles);
      AppHelpers::checkTimeCuts(eventFiles, evtable,
                                m_pars["expcube"], "Exposure");
   }
   m_pars.Prompt("outfile");
   AppHelpers::checkOutputFile(m_pars["clobber"], m_pars["outfile"]);
   m_pars.Prompt("irfs");
   m_pars.Prompt("srcrad");
   m_pars.Prompt("nlong");
   m_pars.Prompt("nlat");
   m_pars.Prompt("nenergies");
   bool compute_submap = m_pars["submap"];
   if (compute_submap) {
      m_pars.Prompt("nlongmin");
      m_pars.Prompt("nlongmax");
      m_pars.Prompt("nlatmin");
      m_pars.Prompt("nlatmax");
   }
   m_pars.Save();
}

void ExpMap::setSourceRegion() {
   m_srRadius = m_pars["srcrad"];
   const RoiCuts & roiCuts = m_helper->observation().roiCuts();
   if (m_srRadius < roiCuts.extractionRegion().radius() + 10.) {
      st_stream::StreamFormatter formatter("gtexpmap", "setSourceRegion", 2);
      formatter.info() << "The radius of the source region, " << m_srRadius 
                       << ", should be significantly larger (say by 10 deg) "
                       << "than the ROI radius of " 
                       << roiCuts.extractionRegion().radius() 
                       << std::endl;
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
   long nlong = m_pars["nlong"];
   long nlat = m_pars["nlat"];
   long nenergies = m_pars["nenergies"];
// Exposure hypercube file.
   std::string expCubeFile = m_pars["expcube"];
   if (expCubeFile != "none") {
      st_facilities::Util::file_ok(expCubeFile);
      ExposureCube & expCube = 
         const_cast<ExposureCube &>(m_helper->observation().expCube());
      expCube.readExposureCube(expCubeFile);
   }
   std::string exposureFile = m_pars["outfile"];
   const Observation & observation = m_helper->observation();
   const RoiCuts & roiCuts = observation.roiCuts();
   bool compute_submap = m_pars["submap"];
   int nlongmin(0);
   int nlongmax(0);
   int nlatmin(0);
   int nlatmax(0);
   if (compute_submap) {
      nlongmin = m_pars["nlongmin"];
      nlongmax = m_pars["nlongmax"];
      nlatmin = m_pars["nlatmin"];
      nlatmax = m_pars["nlatmax"];
   }
   m_helper->observation().expMap().computeMap(exposureFile, observation,
                                               m_srRadius, nlong, nlat,
                                               nenergies, compute_submap,
                                               nlongmin, nlongmax,
                                               nlatmin, nlatmax); 
   std::auto_ptr<tip::Image> 
      image(tip::IFileSvc::instance().editImage(exposureFile, ""));
   roiCuts.writeDssKeywords(image->getHeader());
   roiCuts.writeGtiExtension(exposureFile);
}
