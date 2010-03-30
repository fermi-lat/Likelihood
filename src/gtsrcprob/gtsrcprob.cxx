/**
 * @file gtsrcprob.cxx
 * @brief Add probabilities for each event that it belongs to the
 * various model components given an xml model definition.  These
 * probabilities are proportional to the count rate densities
 * computed in Source::fluxDensity(...).
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/gtsrcprob/gtsrcprob.cxx,v 1.1 2010/03/30 00:03:26 jchiang Exp $
 */

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>

#include "facilities/Util.h"

#include "st_stream/StreamFormatter.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "st_facilities/FitsUtil.h"
#include "st_facilities/Util.h"

#include "xmlBase/Dom.h"
#include "xmlBase/XmlParser.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/DiffRespNames.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/Event.h"
#include "Likelihood/ScData.h"
#include "Likelihood/SourceModel.h"

using namespace Likelihood;

/**
 * @class SourceProbs
 * 
 */

class SourceProbs : public st_app::StApp {

public:

   SourceProbs();

   virtual ~SourceProbs() throw() {
      try {
         delete m_helper;
         delete m_sourceModel;
         delete m_formatter;
      } catch (std::exception & eObj) {
         std::cerr << eObj.what() << std::endl;
      } catch (...) {
         std::cerr << "caught unknown exception in "
                   << "diffuseResponse destructor." 
                   << std::endl;
      }
   }

   virtual void run();
   virtual void banner() const;

private:

   AppHelpers * m_helper;
   SourceModel * m_sourceModel;
   st_stream::StreamFormatter * m_formatter;

   st_app::AppParGroup & m_pars;

   std::vector<std::string> m_srclist;

   void promptForParameters();
   void buildSourceModel();
   void readEventData();
   void writeDensities() const;
   void getSourceList();
   std::string columnName(std::string srcName) const;

   static std::string s_cvs_id;
};

st_app::StAppFactory<SourceProbs> myAppFactory("gtsrcprob");

std::string SourceProbs::s_cvs_id("$Name:  $");

SourceProbs::SourceProbs() 
   : st_app::StApp(), m_helper(0), m_sourceModel(0), 
     m_formatter(new st_stream::StreamFormatter("gtsrcprob", "", 2)),
     m_pars(st_app::StApp::getParGroup("gtsrcprob")) {
   setVersion(s_cvs_id);
}

void SourceProbs::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

void SourceProbs::run() {
   promptForParameters();
   m_helper = new AppHelpers(&m_pars, "UNBINNED");
   m_helper->setRoi("", "EVENTS", false);
   m_helper->readScData();
   buildSourceModel();
   readEventData();
   getSourceList();
   writeDensities();
}

void SourceProbs::promptForParameters() {
   m_pars.Prompt();
   m_pars.Save();
}
 
void SourceProbs::buildSourceModel() {
   m_sourceModel = new SourceModel(m_helper->observation());
   std::string sourceModel = m_pars["srcmdl"];
   st_facilities::Util::file_ok(sourceModel);
   bool requireExposure;
   m_sourceModel->readXml(sourceModel, m_helper->funcFactory(),
                          requireExposure=false);
}

std::string SourceProbs::columnName(std::string srcName) const {
   Event::toLower(srcName);
   return srcName;
}

void SourceProbs::readEventData() {
   std::string evfile = m_pars["evfile"];
   m_helper->observation().eventCont().getEvents(evfile);
}

void SourceProbs::getSourceList() {
   std::string srclist = m_pars["srclist"];
   if (srclist == "" || srclist == "none") {
      m_sourceModel->getSrcNames(m_srclist);
      return;
   }
   std::string skip;
   bool cleanlines;
   st_facilities::Util::readLines(srclist, m_srclist, skip="#",
                                  cleanlines=true);
}

void SourceProbs::writeDensities() const {
   std::string outfile = m_pars["outfile"];

   facilities::Util::expandEnvVar(&outfile);

// Open the output file, copying the FT1 file contents first.
   std::string filter;
   bool clobber;
   st_facilities::FitsUtil::fcopy(m_pars["evfile"], outfile, m_pars["evtable"],
                                  filter="", clobber=m_pars["clobber"]);
   tip::Table * evtable
      = tip::IFileSvc::instance().editTable(outfile, m_pars["evtable"]);

// Add the column names for the probabilities to the output file.
   for (size_t i(0); i < m_srclist.size(); i++) {
      try {
         evtable->appendField(m_srclist.at(i), "1E");
      } catch (tip::TipException & eObj) {
         // Column with this name already exists, so do nothing and 
         // just reuse it.
      }
   }

// Grab the names of all of the model components for the normalization 
// calculation.
   std::vector<std::string> srcNames;
   m_sourceModel->getSrcNames(srcNames);

// Loop over ft1 events and output file rows.
   const std::vector<Event> & 
      events(m_helper->observation().eventCont().events());

   tip::Table::Iterator it = evtable->begin();
   tip::Table::Record & row = *it;
   for (size_t j(0); it != evtable->end(); j++, ++it) {
       std::map<std::string, double> densities;
       double normalization(0);
       // Loop over all sources in the model to determine the normalization.
       for (std::vector<std::string>::iterator name(srcNames.begin());
            name != srcNames.end(); ++name) {
          const Source & src(m_sourceModel->source(*name));
          densities[*name] = src.fluxDensity(events.at(j));
          normalization += densities[*name];
       }
       // Loop over items in the desired source list.
       for (std::vector<std::string>::const_iterator name(m_srclist.begin());
            name != m_srclist.end(); ++name) {
          row[*name].set(densities[*name]/normalization);
       }
   }
   delete evtable;
}
