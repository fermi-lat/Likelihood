/**
 * @file gtwtsmap.cxx
 * @brief Build likelihood weights file.
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

#include "Likelihood/CountsMapBase.h" 
#include "Likelihood/AppHelpers.h" 

/**
 * @class WtsMap
 *
 * @brief Derived class of st_app::StApp for build likelihood weights files.
 */

class WtsMap : public st_app::StApp {

public:

  WtsMap() : st_app::StApp(),
	     m_pars(st_app::StApp::getParGroup("gtwtsmap")),
	     m_alpha_map(0), m_bkg_eff_map(0), m_output(0) {
    setVersion(s_cvs_id);
  }
  virtual ~WtsMap() throw() {
    try {
      delete m_alpha_map;
      delete m_bkg_eff_map;
      delete m_output;
    } catch (std::exception & eObj) {
      std::cerr << eObj.what() << std::endl;
    } catch (...) {
    }
  }
  virtual void run();
  virtual void banner() const;

private:

  st_app::AppParGroup & m_pars;

  Likelihood::CountsMapBase * m_alpha_map;
  Likelihood::CountsMapBase * m_bkg_eff_map;
  Likelihood::CountsMapBase * m_output;

  void computeWtsMap();
  void updateDssKeywords();

  static std::string s_cvs_id;
};

st_app::StAppFactory<WtsMap> myAppFactory("gtwtsmap");

std::string WtsMap::s_cvs_id("$Name:  $");

void WtsMap::banner() const {
  int verbosity = m_pars["chatter"];
  if (verbosity > 2) {
    st_app::StApp::banner();
  }
}

void WtsMap::run() {
  m_pars.Prompt();
  m_pars.Save();
  computeWtsMap();
  updateDssKeywords();
}

void WtsMap::computeWtsMap() {

  std::string effbkgmap_file = m_pars["effbkgfile"];
  std::string alphamap_file = m_pars["alphafile"];
  float epsilon = m_pars["epsilon"];
  float epsilon2 = epsilon*epsilon;

  if ( alphamap_file == "none" ||
       alphamap_file == "" ) {
    m_alpha_map = 0;
  } else {
    m_alpha_map = Likelihood::AppHelpers::readCountsMap(alphamap_file); 
  }
  
  m_bkg_eff_map = Likelihood::AppHelpers::readCountsMap(effbkgmap_file); 
  m_output = Likelihood::CountsMapBase::makeWtsMap(epsilon2, m_alpha_map, *m_bkg_eff_map);
  std::string outfile = m_pars["outfile"];
  m_output->writeAsWeightsMap("gtwtsmap", outfile);
  
}

void WtsMap::updateDssKeywords() {
  Likelihood::CountsMapBase::copyAndUpdateDssKeywords(m_pars["effbkgfile"], 
						      m_pars["outfile"],
						      0,
						      "CALDB");
  Likelihood::CountsMapBase::addWtsMapKeywords(m_pars["outfile"],
					       m_pars["epsilon"],
					       m_pars["effbkgfile"],
					       m_pars["alphafile"]);
}
