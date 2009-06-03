/** 
 * @file gtaddlivetime.cxx
 * @brief Create an Exposure hypercube.
 * @author J. Chiang
 *
 *  $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/gtaddlivetime/gtaddlivetime.cxx,v 1.13 2009/06/02 16:43:44 jchiang Exp $
 */

#include <cstdlib>

#include <iostream>
#include <sstream>
#include <stdexcept>

#include "st_facilities/FitsUtil.h"
#include "st_facilities/Util.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "tip/IFileSvc.h"
#include "tip/Header.h"
#include "tip/Table.h"

#include "dataSubselector/Cuts.h"

#include "Likelihood/AppHelpers.h"

/**
 * @class AddLivetime
 * @brief For two exposure hypercube files, add the livetimes and merge the
 * GTIs.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/gtaddlivetime/gtaddlivetime.cxx,v 1.13 2009/06/02 16:43:44 jchiang Exp $
 */

class AddLivetime : public st_app::StApp {
public:
   AddLivetime() : st_app::StApp(), 
                   m_pars(st_app::StApp::getParGroup("gtltsum")) {
   setVersion(s_cvs_id);
}
   virtual ~AddLivetime() throw() {
      try {
      } catch (std::exception &eObj) {
         std::cerr << eObj.what() << std::endl;
      } catch (...) {
      }
    }
   virtual void run();
   virtual void banner() const;
private:
   st_app::AppParGroup & m_pars;
   std::vector<std::string> m_fileList;

   void promptForParameters();
   void addTables(const std::string & tableName, bool writeGtis=true);
   void writeDateKeywords(const std::string & outfile,
                          double tstart, double tstop) const;

   static std::string s_cvs_id;
};

st_app::StAppFactory<AddLivetime> myAppFactory("gtltsum");

std::string AddLivetime::s_cvs_id("$Name:  $");

void AddLivetime::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

void AddLivetime::run() {
   promptForParameters();
   st_facilities::FitsUtil::fcopy(m_fileList.front(), m_pars["outfile"],
                                  m_pars["table"], "", m_pars["clobber"]);
   addTables(m_pars["table"]);
   try {
      addTables(m_pars["table2"], false);
   } catch(tip::TipException & eObj) {
      std::cout << eObj.what() << std::endl;
   }
}

void AddLivetime::promptForParameters() {
   m_pars.Prompt("infile1");
   std::string infile1 = m_pars["infile1"];
   st_facilities::Util::resolve_fits_files(infile1, m_fileList);
   if (m_fileList.size() == 1) {
      m_pars.Prompt("infile2");
      std::string infile2 = m_pars["infile2"];
      m_fileList.push_back(infile2);
   }
   m_pars.Prompt("outfile");
   m_pars.Save();
}

void AddLivetime::addTables(const std::string & tableName,
                            bool writeGtis) {
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());

   const tip::Table * table(fileSvc.readTable(m_fileList.front(), tableName));
      
   std::vector< std::vector<double> > rows(table->getNumRecords());
   tip::Table::ConstIterator it = table->begin();
   tip::ConstTableRecord & row = *it;
   for (size_t i=0; it != table->end(); ++it, i++) {
      row["COSBINS"].get(rows.at(i));
   }
   const tip::Header & header(table->getHeader());
   double tstart, tstop;
   header["TSTART"].get(tstart);
   header["TSTOP"].get(tstop);
   delete table;

   std::vector<dataSubselector::Cuts> my_cuts;
   my_cuts.push_back(dataSubselector::Cuts(m_fileList.front(), tableName,
                                           false, true));

   for (size_t k=1; k < m_fileList.size(); k++) {
      const tip::Table * table(fileSvc.readTable(m_fileList.at(k), tableName));
      if (rows.size() != static_cast<size_t>(table->getNumRecords())) {
         std::ostringstream message;
         message << "The size of the Exposure extension in " 
                 << m_fileList.at(k)
                 << " does not match the size of the extension in  " 
                 << m_fileList.front();
         throw std::runtime_error(message.str());
      }
      tip::Table::ConstIterator it = table->begin();
      tip::ConstTableRecord & row = *it;
      for (size_t i=0; it != table->end(); ++it, i++) {
         std::vector<double> my_row;
         row["COSBINS"].get(my_row);
         if (rows.at(i).size() != my_row.size()) {
            std::ostringstream message;
            message << "The number of COSBINS columns in "
                    << m_fileList.at(k)
                    << " does not match the number of columns in "
                    << m_fileList.front();
            throw std::runtime_error(message.str());
         }
         for (size_t j=0; j < my_row.size(); j++) {
            rows.at(i).at(j) += my_row.at(j);
         }
      }
      double tvalue;
      table->getHeader()["TSTART"].get(tvalue);
      if (tvalue < tstart) {
         tstart = tvalue;
      }
      table->getHeader()["TSTOP"].get(tvalue);
      if (tvalue > tstop) {
         tstop = tvalue;
      }
      delete table;
      my_cuts.push_back(dataSubselector::Cuts(m_fileList.at(k), tableName,
                                              false, true));
   }

   std::string outfile = m_pars["outfile"];
   tip::Table * outtable(fileSvc.editTable(outfile, tableName));
   tip::Table::Iterator it2 = outtable->begin();
   for (size_t i=0; it2 != outtable->end(); ++it2, i++) {
      (*it2)["COSBINS"].set(rows.at(i));
   }

   dataSubselector::Cuts new_cuts = dataSubselector::Cuts::mergeGtis(my_cuts);
   new_cuts.writeDssKeywords(outtable->getHeader());

   outtable->getHeader()["TSTART"].set(tstart);
   outtable->getHeader()["TSTOP"].set(tstop);

   delete outtable;

   if (writeGtis) {
      new_cuts.writeGtiExtension(outfile);
      writeDateKeywords(outfile, tstart, tstop);
      st_facilities::FitsUtil::writeFilename(outfile);
   }
}

void AddLivetime::writeDateKeywords(const std::string & outfile,
                                    double tstart, double tstop) const {
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   std::vector<std::string> extnames;
   extnames.push_back("");
   extnames.push_back(m_pars["table"]);
   extnames.push_back(m_pars["table2"]);
   extnames.push_back("CTHETABOUNDS");
   extnames.push_back("GTI");
   for (std::vector<std::string>::const_iterator name(extnames.begin());
        name != extnames.end(); ++name) {
      try {
         tip::Extension * hdu(fileSvc.editExtension(outfile, *name));
         st_facilities::Util::writeDateKeywords(hdu, tstart, tstop, *name!="");
         if (*name == "") {
            hdu->getHeader()["CREATOR"].set("gtltsum " + getVersion());
         }
         delete hdu;
      } catch (tip::TipException & eObj) {
         std::cout << eObj.what() << std::endl;
      }
   }
}
