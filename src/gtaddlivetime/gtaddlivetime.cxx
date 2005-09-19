/** 
 * @file gtaddlivetime.cxx
 * @brief Create an Exposure hypercube.
 * @author J. Chiang
 *
 *  $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/gtaddlivetime/gtaddlivetime.cxx,v 1.34 2005/09/12 22:16:34 jchiang Exp $
 */

#include <cstdlib>

#include <iostream>
#include <sstream>
#include <stdexcept>

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "tip/IFileSvc.h"
#include "tip/Header.h"
#include "tip/Table.h"

#include "dataSubselector/Cuts.h"

#include "Likelihood/AppHelpers.h"

#include "Verbosity.h"

/**
 * @class AddLivetime
 * @brief For two exposure hypercube files, add the livetimes and merge the
 * GTIs.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/gtaddlivetime/gtaddlivetime.cxx,v 1.34 2005/09/12 22:16:34 jchiang Exp $
 */

class AddLivetime : public st_app::StApp {
public:
   AddLivetime() : st_app::StApp(), 
                   m_pars(st_app::StApp::getParGroup("gtaddlivetime")) {}
   virtual ~AddLivetime() throw() {
      try {
      } catch (std::exception &eObj) {
         std::cerr << eObj.what() << std::endl;
      } catch (...) {
      }
    }
   virtual void run();
   virtual void banner() const {}
private:
   st_app::AppParGroup & m_pars;
   void promptForParameters();
   void checkGtis();
   void addFiles();
};

st_app::StAppFactory<AddLivetime> myAppFactory("gtaddlivetime");

void AddLivetime::run() {
   promptForParameters();
   checkGtis();
   addFiles();
}

void AddLivetime::promptForParameters() {
   m_pars.Prompt("infile1");
   m_pars.Prompt("infile2");
   m_pars.Prompt("outfile");
   m_pars.Save();
}

void AddLivetime::checkGtis() {
   std::string infile1 = m_pars["infile1"];
   std::string infile2 = m_pars["infile2"];
   std::string table_name = m_pars["table"];
   Likelihood::AppHelpers::checkTimeCuts(infile1, table_name,
                                         infile2, table_name);
}

void AddLivetime::addFiles() {
   std::string infile1 = m_pars["infile1"];
   std::string infile2 = m_pars["infile2"];
   std::string table_name = m_pars["table"];
   const tip::Table * table1 = 
      tip::IFileSvc::instance().readTable(infile1, table_name);
   const tip::Table * table2 = 
      tip::IFileSvc::instance().readTable(infile2, table_name);

   if (table1->getNumRecords() != table2->getNumRecords()) {
      std::ostringstream message;
      message << "The size of the Exposure extension in " << infile1
              << " does not match the size of the extension in  " << infile2;
      throw std::runtime_error(message.str());
   }

   std::string outfile = m_pars["outfile"];
   tip::IFileSvc::instance().createFile(outfile, infile1);
   tip::Table * outtable = 
      tip::IFileSvc::instance().editTable(outfile, table_name);

   tip::Table::ConstIterator it1 = table1->begin();
   tip::ConstTableRecord & row1 = *it1;

   tip::Table::ConstIterator it2 = table2->begin();
   tip::ConstTableRecord & row2 = *it2;

   tip::Table::Iterator it3 = outtable->begin();
   tip::Table::Record & outrow = *it3;

   for ( ; it1 != table1->end(); ++it1, ++it2, ++it3) {
      std::vector<double> x1, x2, sum;
      row1["COSBINS"].get(x1);
      row2["COSBINS"].get(x2);
      for (unsigned int i = 0; i < x1.size(); i++) {
         sum.push_back(x1.at(i) + x2.at(i));
      }
      outrow["COSBINS"].set(sum);
   }

   delete table1;
   delete table2;

   std::vector<dataSubselector::Cuts> my_cuts;

   my_cuts.push_back(dataSubselector::Cuts(infile1, table_name, false, true));
   my_cuts.push_back(dataSubselector::Cuts(infile2, table_name, false, true));

   dataSubselector::Cuts new_cuts =
      dataSubselector::Cuts::mergeGtis(my_cuts);

   tip::Header & my_header(outtable->getHeader());

   double ndskeys;
   my_header["NDSKEYS"].get(ndskeys);
   dataSubselector::Cuts::removeDssKeywords(outfile, table_name, 
                                            static_cast<int>(ndskeys));

   delete outtable;
   outtable = tip::IFileSvc::instance().editTable(outfile, table_name);

   my_header = outtable->getHeader();

   new_cuts.writeDssKeywords(my_header);
   delete outtable;

   new_cuts.writeGtiExtension(outfile);
}
