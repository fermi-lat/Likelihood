/**
 * @file gteffbkg.cxx
 * @brief Compute effective background maps.
 * @author E. Charles
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
#include "Likelihood/CountsMapBase.h" 

/**
 * @class BkgEffMap
 *
 * @brief Derived class of st_app::StApp for making effective background maps.
 * 
 *
 */

class BkgEffMap : public st_app::StApp {

public:

   BkgEffMap() : st_app::StApp(),
                m_pars(st_app::StApp::getParGroup("gteffbkg")),
		 m_helper(0), m_bkg_eff_map(0), m_dataMap(0) {
      setVersion(s_cvs_id);
   }
   virtual ~BkgEffMap() throw() {
      try {
         delete m_helper;
         //delete m_bkg_eff_map;
         //delete m_dataMap;
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
   Likelihood::CountsMapBase * m_dataMap;  
   Likelihood::CountsMapBase * m_bkg_eff_map; 

   void computeBkgEffMap();
   void updateDssKeywords();

   static std::string s_cvs_id;
};

st_app::StAppFactory<BkgEffMap> myAppFactory("gteffbkg");

std::string BkgEffMap::s_cvs_id("$Name:  $");

void BkgEffMap::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

void BkgEffMap::run() {
   m_pars.Prompt();
   m_pars.Save();
   computeBkgEffMap();
   updateDssKeywords();
}

void BkgEffMap::computeBkgEffMap() {

   m_helper = new Likelihood::AppHelpers(&m_pars, "BINNED");
   const Likelihood::MeanPsf& mean_psf = m_helper->observation().meanpsf();
   m_helper->setRoi(m_pars["cmap"], "", false);
   std::string cmapfile = m_pars["cmap"];
   m_dataMap = Likelihood::AppHelpers::readCountsMap(cmapfile);  // EAC: use AppHelpers to read the right type of map
   m_bkg_eff_map = m_dataMap->makeBkgEffMap(mean_psf);
   std::string outfile = m_pars["outfile"];
   m_bkg_eff_map->writeOutput("gteffbkg", outfile);
}

void BkgEffMap::updateDssKeywords() {
  Likelihood::CountsMapBase::copyAndUpdateDssKeywords(m_pars["cmap"], 
						      m_pars["outfile"],
						      m_helper, 
						      m_pars["irfs"]);
}
