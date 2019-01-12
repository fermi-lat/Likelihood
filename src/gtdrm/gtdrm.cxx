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
#include "Likelihood/Drm.h"
#include "Likelihood/FileUtils.h"

/**
 * @class GtDrm
 *
 * @brief Derived class of st_app::StApp for summing up source maps
 * with the spectal fit parameters from a binned likelihood analysis
 * applied.
 *
 */

class GtDrm : public st_app::StApp {

public:

   GtDrm() : st_app::StApp(),
                m_pars(st_app::StApp::getParGroup("gtdrm")),
                m_helper(0), m_logLike(0) {
      setVersion(s_cvs_id);
   }
   virtual ~GtDrm() throw() {
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

   void computeDrm();
   void updateDssKeywords();

   static std::string s_cvs_id;
};

st_app::StAppFactory<GtDrm> myAppFactory("gtdrm");

std::string GtDrm::s_cvs_id("$Name:  $");

void GtDrm::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

void GtDrm::run() {
   m_pars.Prompt();
   m_pars.Save();
   computeDrm();
   updateDssKeywords();
}

void GtDrm::computeDrm() {

   m_helper = new Likelihood::AppHelpers(&m_pars, "BINNED");
   
   m_logLike = Likelihood::AppHelpers::makeBinnedLikelihood(m_pars, *m_helper, "cmap");
   m_dataMap = const_cast<Likelihood::CountsMapBase*>(&(m_logLike->countsMap()));

   const Likelihood::Drm& drm = m_logLike->drm();   

   std::string outfile = m_pars["outfile"];
   m_dataMap->writeEmptyOutput("gtdrm", outfile);
   tip::Extension* table_out = Likelihood::FileUtils::write_drm_to_table(outfile, "DRM", drm);
   delete table_out;
   
}

void GtDrm::updateDssKeywords() {
   Likelihood::CountsMapBase::copyAndUpdateDssKeywords(m_pars["cmap"], 
						       m_pars["outfile"],
						       m_helper, 
						       m_pars["irfs"]);
}
