/**
 * @file gtsrcprob.cxx
 * @brief Add probabilities for each event that it belongs to the
 * various model components given an xml model definition.  These
 * probabilities are proportional to the count rate densities
 * computed in Source::fluxDensity(...).
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/gtsrcprob/gtsrcprob.cxx,v 1.3 2010/04/05 19:13:03 jchiang Exp $
 */

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>

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

   std::string evfile = m_pars["evfile"];
   std::string outfile = m_pars["outfile"];

   if (outfile == evfile) {
      m_formatter->info() << "The output file cannot be the same as the "
                          << "input file. \nPlease specify a different output "
                          << "filename." << std::endl;
      std::exit(0);
   }
   bool clobber = m_pars["clobber"];
   if (!clobber && st_facilities::Util::fileExists(outfile)) {
      m_formatter->info() << "The output file already exists and clobber=yes.\n"
                          << "Please specify a different output "
                          << "filename." << std::endl;
      std::exit(0);
   }      
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

   tip::Table::FieldCont columns(evtable->getValidFields());
// Add the column names for the probabilities to the output file.
   for (size_t i(0); i < m_srclist.size(); i++) {
      // Check if desired column name exists, if not, then add it to
      // the table.
      std::string fieldName(columnName(m_srclist.at(i)));
      if (std::count(columns.begin(), columns.end(), fieldName) == 0) {
         evtable->appendField(m_srclist.at(i), "1E");
         // Delete TNULL that is added incorrectly by tip for floats.
         int fieldIndex = evtable->getFieldIndex(m_srclist.at(i)) + 1;
         std::ostringstream nullkeyword;
         nullkeyword << "TNULL" << fieldIndex;
         try {
            evtable->getHeader().erase(nullkeyword.str());
         } catch (...) {
            // do nothing if tip fails us again
         }
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
