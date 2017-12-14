/**
 * @file gtalphabkg.cxx
 * @brief Compute the factor needed to combine effective background maps.
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

#include "st_facilities/Util.h"

#include "tip/IFileSvc.h"
#include "tip/Image.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/CountsMapBase.h"  

/**
 * @class AlphaBkg
 *
 * @brief Derived class of st_app::StApp for summing up source maps
 * with the spectal fit parameters from a binned likelihood analysis
 * applied.
 *
 */

class AlphaBkg : public st_app::StApp {

public:

   AlphaBkg() : st_app::StApp(),
                m_pars(st_app::StApp::getParGroup("gtalphabkg")),
		m_output(0) {
      setVersion(s_cvs_id);
   }
   virtual ~AlphaBkg() throw() {
      try {
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

   std::vector<Likelihood::CountsMapBase*> m_inputMaps;
   Likelihood::CountsMapBase* m_output;

   void computeAlphaBkg();
   void updateDssKeywords();

   static std::string s_cvs_id;
};

st_app::StAppFactory<AlphaBkg> myAppFactory("gtalphabkg");

std::string AlphaBkg::s_cvs_id("$Name:  $");

void AlphaBkg::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

void AlphaBkg::run() {
   m_pars.Prompt();
   m_pars.Save();
   computeAlphaBkg();
   updateDssKeywords();
}

void AlphaBkg::computeAlphaBkg() {
   
   float epsilon = m_pars["epsilon"];
   float epsilon2 = epsilon*epsilon;

   std::vector<std::string> inputFiles;
   st_facilities::Util::readLines(m_pars["inputs"], inputFiles);

   for ( std::vector<std::string>::const_iterator itr = inputFiles.begin();
	 itr != inputFiles.end(); itr++ ) {
     Likelihood::CountsMapBase* cmap = Likelihood::AppHelpers::readCountsMap(*itr); 
     m_inputMaps.push_back(cmap);
   }

   m_output = Likelihood::CountsMapBase::makeAlphaMap(epsilon2, m_inputMaps);
   std::string outfile = m_pars["outfile"];
   m_output->writeOutput("gtalphgbkg", outfile);

}


void AlphaBkg::updateDssKeywords() {
  std::string outfile = m_pars["outfile"];
  tip::Image * my_image = tip::IFileSvc::instance().editImage(outfile, "");
  my_image->getHeader().setKeyword("NDSKEYS", 0);
}
