/**
 * @file diffuseResponses.cxx
 * @brief Adds diffuse response information for desired components.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/diffuseResponses/diffuseResponses.cxx,v 1.56 2009/02/18 18:13:39 jchiang Exp $
 */

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <fstream>
#include <iostream>

#include "facilities/Util.h"

#include "st_stream/StreamFormatter.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

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

using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;
using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;
using namespace Likelihood;

namespace {
   std::string alt_diffuseSrcName(const std::string srcName,
                                  const Observation & obs) {
      std::string name(obs.respFuncs().respName() + "::" + srcName);
      Event::toLower(name);
      return name;
   }
} // anonymous namespace

/**
 * @class diffuseResponses
 * @brief FTOOL to add diffuse response information to an FT1 file for
 * extragalactic and Galactic diffuse emission.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/diffuseResponses/diffuseResponses.cxx,v 1.56 2009/02/18 18:13:39 jchiang Exp $
 */

class diffuseResponses : public st_app::StApp {

public:

   diffuseResponses();

   virtual ~diffuseResponses() throw() {
      try {
         delete m_helper;
         delete m_srcModel;
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
   SourceModel * m_srcModel;
   st_stream::StreamFormatter * m_formatter;
   
   double m_srRadius;
   st_app::AppParGroup & m_pars;

   std::vector<Event> m_events;
   std::vector<DiffuseSource *> m_srcs;
   std::vector<std::string> m_srcNames;

   DiffRespNames m_columnNames;
   size_t m_ndifrsp;

   bool m_useEdisp;

   void promptForParameters();
   void checkColumnVersion(const std::string & evfile) const;
   void convert_header(const std::string & evfile) const;
   void readDiffuseNames(std::vector<std::string> & srcNames);
   void readDiffRespNames(std::auto_ptr<const tip::Table> events,
                          std::vector<std::string> & colnames);
   void readExistingDiffRespKeys(const tip::Table * events);
   bool haveDiffuseColumns(const std::string & eventFile);
   void buildSourceModel();
   void readEventData(std::string eventFile);
   void computeEventResponses();
   void writeEventResponses(std::string eventFile);
   void getDiffuseSources();
   void setGaussianParams(const Event & event, const std::string & name,
                          tip::Table::Vector<double> & params);
   std::string diffuseSrcName(const std::string & srcName) const;

   static std::string s_cvs_id;
};

st_app::StAppFactory<diffuseResponses> myAppFactory("gtdiffrsp");

std::string diffuseResponses::s_cvs_id("$Name:  $");

diffuseResponses::diffuseResponses() 
   : st_app::StApp(), m_helper(0), m_srcModel(0), 
     m_formatter(new st_stream::StreamFormatter("gtdiffrsp", "", 2)),
     m_srRadius(30.),
     m_pars(st_app::StApp::getParGroup("gtdiffrsp")),
     m_ndifrsp(0) {
   setVersion(s_cvs_id);
}

void diffuseResponses::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

void diffuseResponses::run() {
   promptForParameters();
   bool clobber = m_pars["clobber"];
   m_helper = new AppHelpers(&m_pars, "UNBINNED");
   m_helper->setRoi("", "EVENTS", false);
   m_helper->readScData();
   m_srcModel = new SourceModel(m_helper->observation(), true);
   m_useEdisp = m_pars["edisp"];
   ResponseFunctions & respFuncs = 
      const_cast<ResponseFunctions &>(m_helper->observation().respFuncs());
   respFuncs.setEdispFlag(m_useEdisp);
   std::vector<std::string> eventFiles;
   st_facilities::Util::resolve_fits_files(m_pars["evfile"], eventFiles);
   std::vector<std::string>::const_iterator evtfile;
   buildSourceModel();
   m_formatter->warn() << "Working on...\n";
   for (evtfile = eventFiles.begin(); evtfile != eventFiles.end(); ++evtfile) {
      checkColumnVersion(*evtfile);
      if (clobber || !haveDiffuseColumns(*evtfile)) {
         m_formatter->warn() << *evtfile;
         readEventData(*evtfile);
         computeEventResponses();
         writeEventResponses(*evtfile);
      } else {
         m_formatter->warn() << "Diffuse columns have already been "
                             << "computed for "
                             << *evtfile << "...skipping it." 
                             << std::endl;
      }
   }
}

void diffuseResponses::promptForParameters() {
   m_pars.Prompt();
   m_pars.Save();
}

void diffuseResponses::checkColumnVersion(const std::string & evfile) const {
   const tip::Table *
      events(tip::IFileSvc::instance().readTable(evfile, m_pars["evtable"]));
   const tip::Header & header(events->getHeader());
   if (header.find("NDIFRSP") == header.end()) {
      delete events;
      bool convert = m_pars["convert"];
      if (!convert) {
         std::ostringstream message;
         message <<"NDIFRSP keyword not found in EVENTS HDU of "
                 << evfile
                 << ", and convert=no.\n"
                 << "gtdiffrsp cannot proceed unless you convert this file.";
         throw std::runtime_error(message.str());
      }
      m_formatter->warn() << "Converting EVENTS header for "
                          << evfile
                          << std::endl;
      convert_header(evfile);
      return;
   }
   delete events;
}

void diffuseResponses::convert_header(const std::string & evfile) const {
   std::auto_ptr<tip::Table> 
      events(tip::IFileSvc::instance().editTable(evfile, m_pars["evtable"]));
   tip::Header & header(events->getHeader());
   const tip::Table::FieldCont & validFields(events->getValidFields());
   std::vector<size_t> diffrsp_indx;
   for (size_t i(0); i < validFields.size(); i++) {
      if (validFields.at(i).find("__") != std::string::npos ||
          validFields.at(i).find("::") != std::string::npos) {
         diffrsp_indx.push_back(i);
      }
   }
   header.append(tip::KeyRecord("NDIFRSP", diffrsp_indx.size(), 
                                "Number of diffuse response columns."));
   for (size_t j(0); j < diffrsp_indx.size(); j++) {
      std::ostringstream key;
      key << "DIFRSP" << j;
      header.append(
         tip::KeyRecord(key.str(), 
                        validFields.at(diffrsp_indx.at(j)), 
                        "Diffuse response label for source model component"));
      std::ostringstream ttype;
      ttype << "TTYPE" << diffrsp_indx.at(j) + 1;
      header[ttype.str()].set(key.str());
   }
}

void diffuseResponses::readDiffuseNames(std::vector<std::string> & srcNames) {
   srcNames.clear();
   std::string xmlFile = m_pars["srcmdl"];
   xmlBase::XmlParser * parser = new xmlBase::XmlParser();
   DOMDocument * doc = parser->parse(xmlFile.c_str());
   DOMElement * source_library = doc->getDocumentElement();
   std::vector<DOMElement *> srcs;
   xmlBase::Dom::getChildrenByTagName(source_library, "source", srcs);
   for (unsigned int i = 0; i < srcs.size(); i++) {
      if (xmlBase::Dom::getAttribute(srcs[i], "type") == "DiffuseSource") {
         srcNames.push_back(xmlBase::Dom::getAttribute(srcs[i], "name"));
      }
   }
}

bool diffuseResponses::haveDiffuseColumns(const std::string & eventFile) {
   std::auto_ptr<const tip::Table> 
      events(tip::IFileSvc::instance().readTable(eventFile,
                                                 m_pars["evtable"]));
   std::vector<std::string> colNames;
   readDiffRespNames(events, colNames);

   std::vector<std::string> srcNames;
   readDiffuseNames(srcNames);
   for (std::vector<std::string>::iterator name = srcNames.begin();
        name != srcNames.end(); ++name) {
      std::string altName(::alt_diffuseSrcName(*name, m_helper->observation()));
      *name = diffuseSrcName(*name);
      if (std::find(colNames.begin(), colNames.end(), *name) == colNames.end()
          && std::find(colNames.begin(), colNames.end(), altName) 
          == colNames.end()) {
         return false;
      }
   }
   return true;
}

void diffuseResponses::
readDiffRespNames(std::auto_ptr<const tip::Table> events,
                  std::vector<std::string> & colnames) {
// Read in DIFRSPxx names (or column names if using old format)
   const tip::Header & header(events->getHeader());
   int nkeys;
   try {
      header["NDIFRSP"].get(nkeys);
      m_ndifrsp = nkeys;
   } catch(tip::TipException & eObj) {
// Assume we have the old format, so just read in column names.
      colnames = events->getValidFields();
      return;
   }
   for (size_t i(0); i < m_ndifrsp; i++) {
      std::ostringstream keyname;
      keyname << "DIFRSP" << i;
      std::string colname;
      header[keyname.str()].get(colname);
      colnames.push_back(colname);
   }
}

std::string diffuseResponses::
diffuseSrcName(const std::string & srcName) const {
   std::string name(m_helper->observation().respFuncs().respName() +
                    "__" + srcName);
   Event::toLower(name);
   return name;
}

void diffuseResponses::buildSourceModel() {
   std::string sourceModel = m_pars["srcmdl"];
   st_facilities::Util::file_ok(sourceModel);
   m_srcModel->readXml(sourceModel, m_helper->funcFactory(), false, false);
}

void diffuseResponses::readEventData(std::string eventFile) {
   m_events.clear();
   facilities::Util::expandEnvVar(&eventFile);
   st_facilities::Util::file_ok(eventFile);
   const tip::Table * events 
      = tip::IFileSvc::instance().readTable(eventFile, m_pars["evtable"]);

   int evclsver(0); // version of event class definition

   const tip::Header & header(events->getHeader());
   try {
      header["EVCLSVER"].get(evclsver);
   } catch(tip::TipException) {
      // keyword missing so use default value
   }

   double ra;
   double dec;
   double energy;
   double time;
   double zenAngle;
   int event_class;
   int conversion_type;
   int eventType;

   ScData & scData = const_cast<ScData &>(m_helper->observation().scData());

   tip::Table::ConstIterator it = events->begin();
   tip::ConstTableRecord & event = *it;
   for ( ; it != events->end(); ++it) {
      event["ra"].get(ra);
      event["dec"].get(dec);
      event["energy"].get(energy);
      event["time"].get(time);
      event["zenith_angle"].get(zenAngle);
      event["event_class"].get(event_class);
      event["conversion_type"].get(conversion_type);
      if (evclsver == 0) {
         eventType = conversion_type;
      } else {
         eventType = conversion_type + 2*event_class;
      }
      Event thisEvent(ra, dec, energy, time, scData.zAxis(time), 
                      scData.xAxis(time), cos(zenAngle*M_PI/180.), 
                      m_helper->observation().respFuncs().useEdisp(),
                      m_helper->observation().respFuncs().respName(),
                      eventType);
      thisEvent.set_classLevel(event_class);
      m_events.push_back(thisEvent);
   }
   delete events;
}

void diffuseResponses::computeEventResponses() {
   getDiffuseSources();
   std::vector<Event>::iterator it = m_events.begin();
   int classLevel_min = m_pars["evclsmin"];
   for (int i = 0; it != m_events.end(); ++it, i++) {
      int factor(m_events.size()/20);
      if (factor == 0) {
         factor = 1;
      }
      if ((i % factor) == 0) {
         m_formatter->warn() << ".";
      }
/// @todo Implement an accurate, faster default calculation; use Gaussian
/// quadrature version for now.
//       it->computeResponse(m_srcs, m_helper->observation().respFuncs(), 
//                           m_srRadius);
      bool useDummyValue(it->classLevel() < classLevel_min);
      it->computeResponseGQ(m_srcs, m_helper->observation().respFuncs(),
                            useDummyValue);
   }
   m_formatter->warn() << "!" << std::endl;
}

void diffuseResponses::writeEventResponses(std::string eventFile) {
   if (m_srcNames.size() == 0) {
      return;
   }
   facilities::Util::expandEnvVar(&eventFile);
   st_facilities::Util::file_ok(eventFile);
   tip::Table * events 
      = tip::IFileSvc::instance().editTable(eventFile, m_pars["evtable"]);
   if (static_cast<size_t>(events->getNumRecords()) != m_events.size()) {
      throw std::runtime_error("diffuseResponses::writeEventResponses:" + 
                               ("\nNumber of records in " 
                                + eventFile 
                                + " does not match number of events."));
   }
   readExistingDiffRespKeys(events);
/// Add the column names to m_columnNames.
   for (size_t i(0); i < m_srcNames.size(); i++) {
      try {
         std::string columnName(diffuseSrcName(m_srcNames[i]));
         m_columnNames.addColumn(columnName);
         std::string fieldName(m_columnNames.key(columnName));
         if (m_useEdisp) {
// Add a 3 dim vector containing the Gaussian parameters describing
// the energy response.
            events->appendField(fieldName, "3E");
         } else {
// Infinite energy response, so just add the single value.
            events->appendField(fieldName, "1E");
         }
// Repair field by removing incorrect TNULL keyword that is added by tip:
         int fieldIndex = events->getFieldIndex(fieldName) + 1;
         std::ostringstream nullkeyword;
         nullkeyword << "TNULL" << fieldIndex;
         try {
            events->getHeader().erase(nullkeyword.str());
         } catch (...) {
            // do nothing if tip fails us again here.
         }
      } catch (tip::TipException &eObj) {
//          m_formatter->warn() << eObj.what() << "\n"
//                              << "Using existing column." << std::endl;
      }
   }
   tip::Table::Iterator it = events->begin();
   tip::Table::Record & row = *it;
   for (int j = 0 ; it != events->end(); j++, ++it) {
      std::vector<std::string>::iterator name = m_srcNames.begin();
      for ( ; name != m_srcNames.end(); ++name) {
         std::string fieldName = m_columnNames.key(diffuseSrcName(*name));
         if (m_useEdisp) {
            tip::Table::Vector<double> respParams = row[fieldName];
            setGaussianParams(m_events[j], *name, respParams);
         } else {
// Assume infinite energy resolution.
            row[fieldName].set(m_events[j].diffuseResponse(1., *name));
         }
      }
   }
// Set DIFRSPxx keywords.
   for (size_t i(0); i < m_columnNames.size(); i++) {
      std::string diffrspName(m_columnNames[i]);
      std::string key(m_columnNames.key(diffrspName));
      events->getHeader()[key].set(diffrspName);
   }
   if (m_columnNames.size() > m_ndifrsp) {
      events->getHeader()["NDIFRSP"].set(m_columnNames.size());
   }      
   delete events;
}

void diffuseResponses::setGaussianParams(const Event & event,
                                         const std::string & name,
                                         tip::Table::Vector<double> & params) {
   double norm, mean, sigma;
   event.computeGaussianParams(name, norm, mean, sigma);
   params[0] = norm;
   params[1] = mean;
   params[2] = sigma;
}

void diffuseResponses::getDiffuseSources() {
   m_srcs.clear();
   m_srcNames.clear();
   std::vector<std::string> srcNames;
   m_srcModel->getSrcNames(srcNames);
   for (unsigned int i = 0; i < srcNames.size(); i++) {
      Source * my_src = m_srcModel->getSource(srcNames[i]);
      if (my_src->getType() == std::string("Diffuse")) {
         m_srcs.push_back(dynamic_cast<DiffuseSource *>(my_src));
         m_srcNames.push_back(srcNames[i]);
      }
   }
}

void diffuseResponses::readExistingDiffRespKeys(const tip::Table * events) {
   const tip::Header & header(events->getHeader());
   int nkeys;
   header["NDIFRSP"].get(nkeys);
   for (int i(0); i < nkeys; i++) {
      std::ostringstream keyname;
      keyname << "DIFRSP" << i;
      std::string colname;
      header[keyname.str()].get(colname);
      if (colname != "NONE") {
         m_columnNames.addColumn(colname);
      }
   }
}
