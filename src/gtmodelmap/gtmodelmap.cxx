/**
 * @file gtmodelmap.cxx
 * @brief Compute a model counts map based on binned likelihood fits.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/gtmodelmap/gtmodelmap.cxx,v 1.46 2016/08/05 21:04:44 echarles Exp $
 */

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "st_stream/StreamFormatter.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "tip/IFileSvc.h"
#include "tip/Image.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/BinnedLikelihood.h"
#include "Likelihood/CountsMapBase.h"  // EAC: switch to new base class
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
   Likelihood::CountsMapBase * m_dataMap;  // EAC: switch to new base class
   Likelihood::BinnedLikelihood * m_logLike;

   void computeModelMap();
   void updateDssKeywords();

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
   computeModelMap();
   updateDssKeywords();
}

void ModelMap::computeModelMap() {

   m_helper = new Likelihood::AppHelpers(&m_pars, "BINNED");
   
   m_logLike = Likelihood::AppHelpers::makeBinnedLikelihood(m_pars, *m_helper, "srcmaps");
   m_dataMap = const_cast<Likelihood::CountsMapBase*>(&(m_logLike->countsMap()));

   m_logLike->set_use_single_fixed_map(false);   
   bool requireExposure, addPointSources, loadMaps, createAllMaps;
   m_logLike->readXml(m_pars["srcmdl"], m_helper->funcFactory(),
                      requireExposure=false, addPointSources=true,
                      loadMaps=false);

   std::vector<float> ext_model_map;
   m_logLike->computeModelMap(ext_model_map);

   Likelihood::ModelMap modelMap(*m_logLike, &ext_model_map);
   
   std::string outfile = m_pars["outfile"];
   std::string outtype = m_pars["outtype"];

   modelMap.writeOutputMap(outfile, outtype);
}

void ModelMap::updateDssKeywords() {
   Likelihood::CountsMapBase::copyAndUpdateDssKeywords(m_pars["srcmaps"], 
						       m_pars["outfile"],
						       m_helper, 
						       m_pars["irfs"]);
   std::string phased_expmap = m_pars["phased_expmap"];
   if (phased_expmap != "none" && phased_expmap != "") {
     // EAC, add the phased exposure map filename to the head keywords
     Likelihood::CountsMapBase::addPhasedExpMapKeywords(m_pars["outfile"], m_pars["phased_expmap"]);
   }
}
