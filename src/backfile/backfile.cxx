/**
 * @file backfile.cxx
 * @brief Application to create "background" pha file for an Xspec 
 * point-source analysis using a fitted likelihood model (that excludes
 * the source in question).
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/backfile/backfile.cxx,v 1.12 2009/12/16 19:05:46 elwinter Exp $
 */

#include <cstdlib>

#include <algorithm>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

#include "st_stream/StreamFormatter.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"
#include "tip/TipFile.h"

#include "st_facilities/FitsUtil.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/LogLike.h"

/**
 * @class BackFile
 * @brief Derived class of st_app::StApp used to compute a background 
 * pha file.
 * 
 */

class BackFile : public st_app::StApp {
public:
   BackFile() : st_app::StApp(),
                m_pars(st_app::StApp::getParGroup("gtbkg")) {
      setVersion(s_cvs_id);
   }
   virtual ~BackFile() throw() {}

   virtual void run();
   virtual void banner() const;
private:
   st_app::AppParGroup & m_pars;
   Likelihood::AppHelpers * m_helper;

   void setup();
   void getEbounds(std::vector<double> & emin,
                   std::vector<double> & emax) const;
   void updateKeywords() const;
   double NpredError(const Likelihood::LogLike & logLike,
                     const std::vector<std::string> & source_names,
                     double emin, double emax) const;
   void writeBackFile(const std::vector<double> & bg_rate,
                      const std::vector<double> & bg_rate_err) const;
   void setHeaderKeyword(const std::string & phafile, 
                         const std::string & extension,
                         const std::string & keyname,
                         const std::string & value) const;
   static std::string s_cvs_id;
};

st_app::StAppFactory<BackFile> myAppFactory("gtbkg");

std::string BackFile::s_cvs_id("$Name:  $");

void BackFile::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

void BackFile::setup() {
   m_pars.Prompt();
   m_pars.Save();
   m_helper = new Likelihood::AppHelpers(&m_pars, "none");
   std::string expMap = m_pars["expmap"];
   m_helper->observation().expMap().readExposureFile(expMap);
//   m_helper->observation().roiCuts().readCuts(expMap, "");
   m_helper->setRoi(expMap, "");
   m_helper->readScData();
//    std::string scfile = m_pars["scfile"];
//    m_helper->observation().scData().readData(scfile);
   std::string expCube = m_pars["expcube"];
   m_helper->observation().expCube().readExposureCube(expCube);

   std::string phafile = m_pars["phafile"];
   m_helper->checkCuts(phafile, "SPECTRUM", expMap, "");
}

void BackFile::run() {
   setup();

   std::string infile = m_pars["phafile"];
   std::string outfile = m_pars["outfile"];

   tip::TipFile inputfile = tip::IFileSvc::instance().openFile(infile);

   inputfile.copyFile(outfile);

   Likelihood::LogLike logLike(m_helper->observation());
   std::string srcModel = m_pars["srcmdl"];
   logLike.readXml(srcModel, m_helper->funcFactory());

   std::vector<std::string> srcNames;
   logLike.getSrcNames(srcNames);

   std::string target = m_pars["target"];

   st_stream::StreamFormatter formatter("gtbkg", "run", 2);

   if (std::find(srcNames.begin(), srcNames.end(), target)
       != srcNames.end()) {
      formatter.info() << "Excluding source " << target 
                       << " from background model." << std::endl;
   } else if (target != "none" && target != "") {
      formatter.info() << "Source named '" << target << "' not found.\n"
                       << "Using all sources in input model for "
                       << "background estimate." << std::endl;
   }

   std::vector<double> emin;
   std::vector<double> emax;
   getEbounds(emin, emax);

   const std::vector<double> & energies = 
      m_helper->observation().roiCuts().energies();

   // Read the exposure time to convert counts to rate.
   const tip::Table * table = tip::IFileSvc::instance().readTable(infile,
                                                                  "SPECTRUM");
   double exposure;
   table->getHeader()["EXPOSURE"].get(exposure);
   delete table;

   std::vector<double> bg_rate;
   std::vector<double> bg_rate_err;
   for (unsigned int i = 0; i < emin.size(); i++) {
      double emin_val = std::max(emin.at(i), energies.front());
      double emax_val = std::min(emax.at(i), energies.back());
      double counts(0);
      for (std::vector<std::string>::const_iterator srcName = srcNames.begin();
           srcName != srcNames.end(); ++srcName) {
         if (*srcName != target) {
            counts += logLike.getSource(*srcName)->Npred(emin_val, emax_val);
         }
      }
      bg_rate.push_back(counts/exposure);
      bg_rate_err.push_back(NpredError(logLike, srcNames, emin_val, emax_val)/
                            exposure);
   }

   updateKeywords();
   writeBackFile(bg_rate, bg_rate_err);
   setHeaderKeyword(infile, "SPECTRUM", "BACKFILE", outfile);

   st_facilities::FitsUtil::writeChecksums(outfile);
}

double BackFile::NpredError(const Likelihood::LogLike & logLike,
                            const std::vector<std::string> & source_names,
                            double emin, double emax) const {
   std::string target = m_pars["target"];
   double variance(0);
   for (size_t i(0); i < source_names.size(); i++) {
      if (source_names[i] == target) {
         continue;
      }
      std::vector<std::string> parnames;
      const Likelihood::Source & src(logLike.source(source_names[i]));
      src.spectrum().getFreeParamNames(parnames);
      for (size_t j(0); j < parnames.size(); j++) {
         double dNpred_dpar(src.NpredDeriv(parnames[j], emin, emax));
         double par_error(src.spectrum().parameter(parnames[j]).error());
         variance += dNpred_dpar*dNpred_dpar*par_error*par_error;
      }
   }
   return std::sqrt(variance);
}

void BackFile::getEbounds(std::vector<double> & emin,
                          std::vector<double> & emax) const {
   std::string pha_file = m_pars["phafile"];
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

// Convert COUNTS to RATE so that floating point values may be used
void BackFile::updateKeywords() const {
   tip::Table * my_spectrum = 
      tip::IFileSvc::instance().editTable(m_pars["outfile"], "SPECTRUM");
   tip::Header & header(my_spectrum->getHeader());
   header["TTYPE2"].set("RATE");
   header["TFORM2"].set("E");
   delete my_spectrum;
}

void BackFile::writeBackFile(const std::vector<double> & bg_rate,
                             const std::vector<double> & bg_rate_err) const {
   std::string outfile = m_pars["outfile"];
   tip::Table * my_spectrum = 
      tip::IFileSvc::instance().editTable(outfile, "SPECTRUM");
   tip::Table::Iterator it = my_spectrum->begin();
   tip::Table::Record & row = *it;
   for (unsigned int i = 0; it != my_spectrum->end(); ++it, i++) {
      row["RATE"].set(bg_rate.at(i));
      row["STAT_ERR"].set(bg_rate_err.at(i));
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
