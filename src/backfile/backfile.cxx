/**
 * @file backfile.cxx
 * @brief Application to create "background" pha file for an Xspec 
 * point-source analysis using a fitted likelihood model (that excludes
 * the source in question).
 * @author J. Chiang
 *
 * $Header$
 */

#include <cstdlib>

#include <iostream>
#include <sstream>
#include <stdexcept>

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "src/FitsTipFile.h"
#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/LogLike.h"

/**
 * @class BackFile
 * @brief Derived class of st_app::StApp used to compute a background 
 * pha file.
 * 
 * @author J. Chiang
 *
 * $Header$
 */

class BackFile : public st_app::StApp {
public:
   BackFile() : st_app::StApp(),
                m_pars(st_app::StApp::getParGroup("gtbackfile")) {}
   virtual ~BackFile() throw() {}

   virtual void run();
   virtual void banner() const {}
private:
   st_app::AppParGroup & m_pars;
   Likelihood::AppHelpers * m_helper;
   void setup();
   void getEbounds(std::vector<double> & emin,
                   std::vector<double> & emax) const;
   void writeBackFile(const std::vector<double> & bg_counts) const;
   void setHeaderKeyword(const std::string & phafile, 
                         const std::string & extension,
                         const std::string & keyname,
                         const std::string & value) const;
};

st_app::StAppFactory<BackFile> myAppFactory("gtbackfile");

void BackFile::setup() {
   m_pars.Prompt();
   m_pars.Save();
   m_helper = new Likelihood::AppHelpers(&m_pars);
   std::string scfile = m_pars["scfile"];
   m_helper->observation().scData().readData(scfile);
   std::string expMap = m_pars["exposure_map_file"];
   m_helper->observation().expMap().readExposureFile(expMap);
   m_helper->observation().roiCuts().readCuts(expMap, "");
   std::string expCube = m_pars["exposure_cube_file"];
   m_helper->observation().expCube().readExposureCube(expCube);
}

void BackFile::run() {
   setup();
   Likelihood::LogLike logLike(m_helper->observation());
   std::string srcModel = m_pars["source_model_file"];
   logLike.readXml(srcModel, m_helper->funcFactory());

   std::vector<std::string> srcNames;
   logLike.getSrcNames(srcNames);

   std::string target = m_pars["target_source"];
   unsigned int chatter = m_pars["chatter"];
   if (chatter > 1) {
      if (std::find(srcNames.begin(), srcNames.end(), target)
          != srcNames.end()) {
         std::cout << "Excluding source " << target 
                   << " from background model.";
      } else if (target != "none" && target != "") {
         std::cout << "Source named '" << target << "' not found.\n"
                   << "Using all sources in input model for "
                   << "background estimate.";
      }
   }         

   std::vector<double> emin;
   std::vector<double> emax;
   getEbounds(emin, emax);

   const std::vector<double> & energies = 
      m_helper->observation().roiCuts().energies();

   std::vector<double> bg_counts;
   for (unsigned int i = 0; i < emin.size(); i++) {
      double counts(0);
      for (std::vector<std::string>::const_iterator srcName = srcNames.begin();
           srcName != srcNames.end(); ++srcName) {
         if (*srcName != target) {
            double emin_val = std::max(emin.at(i), energies.front());
            double emax_val = std::min(emax.at(i), energies.back());
            counts += logLike.getSource(*srcName)->Npred(emin_val, emax_val);
         }
      }
      bg_counts.push_back(counts);
   }
   
   std::string infile = m_pars["pha_file"];
   std::string outfile = m_pars["outfile"];

   tip::FitsTipFile inputfile(infile);
   inputfile.copyFile(outfile);

   writeBackFile(bg_counts);
   setHeaderKeyword(infile, "SPECTRUM", "BACKFILE", outfile);
}

void BackFile::getEbounds(std::vector<double> & emin,
                          std::vector<double> & emax) const {
   std::string pha_file = m_pars["pha_file"];
   std::auto_ptr<const tip::Table>
      ebounds(tip::IFileSvc::instance().readTable(pha_file, "EBOUNDS"));

   long nenergies = ebounds->getNumRecords();
   emin.resize(nenergies);
   emax.resize(nenergies);

   tip::Table::ConstIterator it = ebounds->begin();
   tip::ConstTableRecord & row = *it;
   
   for (unsigned int i = 0; it != ebounds->end(); ++it, i++) {
      row["E_MIN"].get(emin.at(i));
      emin.at(i) /= 1e3;
      row["E_MAX"].get(emax.at(i));
      emax.at(i) /= 1e3;
   }
}

void BackFile::writeBackFile(const std::vector<double> & bg_counts) const {
   std::string outfile = m_pars["outfile"];
   tip::Table * my_spectrum = 
      tip::IFileSvc::instance().editTable(outfile, "SPECTRUM");
   tip::Table::Iterator it = my_spectrum->begin();
   tip::Table::Record & row = *it;
   for (unsigned int i = 0; it != my_spectrum->end(); ++it, i++) {
      row["COUNTS"].set(bg_counts.at(i));
   }
   delete my_spectrum;
}

void BackFile::setHeaderKeyword(const std::string & phafile, 
                                const std::string & extension,
                                const std::string & keyname,
                                const std::string & value) const {
   std::auto_ptr<tip::Table>
      spectrum(tip::IFileSvc::instance().editTable(phafile, extension));
   tip::Header & header = spectrum->getHeader();
   header[keyname].set(value);
}
