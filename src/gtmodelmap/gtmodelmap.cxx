/**
 * @file gtmodelmap.cxx
 * @brief Compute a model counts map based on binned likelihood fits.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/gtmodelmap/gtmodelmap.cxx,v 1.36 2012/09/13 19:24:31 jchiang Exp $
 */

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "st_stream/StreamFormatter.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "dataSubselector/Cuts.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/BinnedLikelihood.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/ModelMap.h"

/**
 * @class ModelMap
 *
 * @brief Derived class of st_app::StApp for summing up source maps
 * with the spectal fit parameters from a binned likelihood analysis
 * applied.
 *
 */

class ModelMap : public st_app::StApp {

public:

   ModelMap() : st_app::StApp(),
                m_pars(st_app::StApp::getParGroup("gtmodel")),
                m_helper(0), m_logLike(0) {
      setVersion(s_cvs_id);
   }
   virtual ~ModelMap() throw() {
      try {
         delete m_logLike;
         delete m_dataMap;
         delete m_helper;
      } catch (std::exception & eObj) {
         std::cerr << eObj.what() << std::endl;
      } catch (...) {
      }
   }
   virtual void run();
   virtual void banner() const;

private:

   st_app::AppParGroup & m_pars;

   Likelihood::AppHelpers * m_helper;
   Likelihood::CountsMap * m_dataMap;
   Likelihood::BinnedLikelihood * m_logLike;

   void computeModelMap();

   static std::string s_cvs_id;
};

st_app::StAppFactory<ModelMap> myAppFactory("gtmodel");

std::string ModelMap::s_cvs_id("$Name:  $");

void ModelMap::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

void ModelMap::run() {
   m_pars.Prompt();
   m_pars.Save();

   std::string smaps = m_pars["srcmaps"];
   std::string irfs = m_pars["irfs"];
   dataSubselector::Cuts::checkIrfs(smaps, "", irfs);

   computeModelMap();
}

void ModelMap::computeModelMap() {
   m_helper = new Likelihood::AppHelpers(&m_pars, "BINNED");
   m_helper->observation().expCube().readExposureCube(m_pars["expcube"]);
   m_helper->setRoi(m_pars["srcmaps"], "", false);
   std::string cmapfile = m_pars["srcmaps"];
   m_dataMap = new Likelihood::CountsMap(cmapfile);
   bool computePointSources, apply_psf_corrections;
   bool performConvolution = m_pars["convol"];
   bool resample = m_pars["resample"];
   int resamp_factor = m_pars["rfactor"];
   double rfactor = static_cast<double>(resamp_factor);
   m_logLike = new Likelihood::BinnedLikelihood(*m_dataMap,
                                                m_helper->observation(),
                                                cmapfile, 
                                                computePointSources=true, 
                                                apply_psf_corrections=true,
                                                performConvolution,
                                                resample, rfactor);
   std::string bexpmap = m_pars["bexpmap"];
   Likelihood::AppHelpers::checkExposureMap(cmapfile, bexpmap);
   bool requireExposure, addPointSources, loadMaps, createAllMaps;
   m_logLike->readXml(m_pars["srcmdl"], m_helper->funcFactory(),
                      requireExposure=false, addPointSources=true,
                      loadMaps=false, createAllMaps=true);

   Likelihood::ModelMap modelMap(*m_logLike);
   
   std::string outfile = m_pars["outfile"];
   std::string outtype = m_pars["outtype"];

   modelMap.writeOutputMap(outfile, outtype);
}
