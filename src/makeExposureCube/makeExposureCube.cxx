/** 
 * @file makeExposureCube.cxx
 * @brief Create an Exposure hypercube.
 * @author J. Chiang
 *
 *  $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/makeExposureCube/makeExposureCube.cxx,v 1.60 2012/10/02 22:09:18 jchiang Exp $
 */

#include <cstdio>
#include <cstdlib>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <memory>
#include <stdexcept>

#include "st_stream/StreamFormatter.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "st_facilities/Util.h"

#include "healpix/CosineBinner.h"

#include "Likelihood/LikeExposure.h"
#include "Likelihood/RoiCuts.h"

namespace {
   void getTBounds(const std::vector< std::pair<double, double> > & gtis,
                   double & tmin, double & tmax) {
      if (gtis.size() > 0) {
         tmin = gtis.at(0).first;
         tmax = gtis.at(0).second;
      }
      for (size_t i(0); i < gtis.size(); i++) {
         if (gtis.at(i).first < tmin) {
            tmin = gtis.at(i).first;
         }
         if (gtis.at(i).second > tmax) {
            tmax = gtis.at(i).second;
         }
      }         
   }
}

/**
 * @class ExposureCube
 * @brief Class to encapsulate methods for creating an exposure
 * hypercube in (ra, dec, cos_theta) using the LikeExposure class.
 *
 * @author J. Chiang
 */
class ExposureCube : public st_app::StApp {
public:
   ExposureCube() : st_app::StApp(), 
                    m_pars(st_app::StApp::getParGroup("gtltcube")), 
                    m_exposure(0), m_roiCuts(0) {
      setVersion(s_cvs_id);
   }
   virtual ~ExposureCube() throw() {
      try {
         delete m_exposure;
         delete m_roiCuts;
      } catch (std::exception &eObj) {
         std::cerr << eObj.what() << std::endl;
      } catch (...) {
      }
    }
   virtual void run();
   virtual void banner() const;
private:
   st_app::AppParGroup & m_pars;
   Likelihood::LikeExposure * m_exposure;
   Likelihood::RoiCuts * m_roiCuts;
   void readRoiCuts();
   void createDataCube();
   void writeTableKeywords(const std::string & outfile,
                           const std::string & tablename) const;
   void writeDateKeywords(const std::string & outfile, 
                          double tstart, double tstop) const;
   static std::string s_cvs_id;
};

st_app::StAppFactory<ExposureCube> myAppFactory("gtltcube");

std::string ExposureCube::s_cvs_id("$Name:  $");

void ExposureCube::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

void ExposureCube::run() {
   m_pars.Prompt();
   m_pars.Save();
   readRoiCuts();
   std::string output_file = m_pars["outfile"];
   if (st_facilities::Util::fileExists(output_file)) {
      if (m_pars["clobber"]) {
         std::remove(output_file.c_str());
      } else {
         st_stream::StreamFormatter formatter("gtltcube", "run", 2);
         formatter.err() << "Output file " << output_file 
                         << " already exists,\n"
                         << "and you have set 'clobber' to 'no'.\n"
                         << "Please provide a different output file name." 
                         << std::endl;
         std::exit(1);
      }
   }
   createDataCube();
   m_exposure->writeFile(output_file);
   writeTableKeywords(output_file, "EXPOSURE");
   writeTableKeywords(output_file, "WEIGHTED_EXPOSURE");

   m_roiCuts->writeGtiExtension(output_file);
   double tstart(m_roiCuts->minTime());
   double tstop(m_roiCuts->maxTime());

   writeDateKeywords(output_file, tstart, tstop);
}

void ExposureCube::writeTableKeywords(const std::string & outfile,
                                      const std::string & tablename) const {
   std::unique_ptr<tip::Table> 
      table(tip::IFileSvc::instance().editTable(outfile, tablename));
   m_roiCuts->writeDssTimeKeywords(table->getHeader());
   double tstart(m_roiCuts->minTime());
   double tstop(m_roiCuts->maxTime());
   tip::Header & header(table->getHeader());
   header["TSTART"].set(tstart);
   header["TSTOP"].set(tstop);
   double zmax = m_pars["zmax"];
   double zmin = m_pars["zmin"];
   if (zmax < 180.) {
      header["ZENMAX"].set(zmax);
      header["ZENMIN"].set(zmin);
   }
}

void ExposureCube::writeDateKeywords(const std::string & outfile, 
                                     double tstart, double tstop) const {
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   std::vector<std::string> extnames;
   extnames.push_back("");
   extnames.push_back("EXPOSURE");
   extnames.push_back("WEIGHTED_EXPOSURE");
   extnames.push_back("CTHETABOUNDS");
   extnames.push_back("GTI");
   for (std::vector<std::string>::const_iterator name(extnames.begin());
        name != extnames.end(); ++name) {
      tip::Extension * hdu(fileSvc.editExtension(outfile, *name));
      st_facilities::Util::writeDateKeywords(hdu, tstart, tstop, *name!="");
      if (*name == "") {
         hdu->getHeader()["CREATOR"].set("gtltcube " + getVersion());
         std::string file_version = m_pars["file_version"];
         hdu->getHeader()["VERSION"].set(file_version);
      }
      delete hdu;
   }
}

void ExposureCube::readRoiCuts() {
   std::string event_file = m_pars["evfile"];
   m_roiCuts = new Likelihood::RoiCuts();
   if (event_file != "") {
      std::string evtable = m_pars["evtable"];
      std::vector<std::string> eventFiles;
      st_facilities::Util::resolve_fits_files(event_file, eventFiles);
      m_roiCuts->readCuts(eventFiles, evtable, false);
   } else {
      double tmin = m_pars["tmin"];
      double tmax = m_pars["tmax"];
      m_roiCuts->setCuts(0, 0, 180, 30, 3e5, tmin, tmax, -1, true);
   }
}

void ExposureCube::createDataCube() {
   st_stream::StreamFormatter formatter("gtltcube", 
                                        "createDataCube", 2);

   std::vector<std::pair<double, double> > gtis;
   m_roiCuts->getGtis(gtis);

   double tmin, tmax;
   ::getTBounds(gtis, tmin, tmax);
   static double maxIntervalSize(30);
   tmin -= 2.*maxIntervalSize;
   tmax += 2.*maxIntervalSize;
   std::ostringstream filter;
   filter << std::setprecision(20);
   filter << "(START >= " << tmin << ") && (STOP <= " << tmax << ")";
   formatter.info(4) << "applying filter: " << filter.str() << std::endl;

   double zmax = m_pars["zmax"];
   double zmin = m_pars["zmin"];

   // Set the number of phibins using the static function interface
   // from healpix::CosineBinner (this is how
   // map_tools/exposure_cube.cxx does it.)
   long nphibins = m_pars["phibins"];
   if (nphibins > 0) {
      healpix::CosineBinner::setPhiBins(nphibins);
   }

   // timeCuts is just needed for LikeExposure constructor, but it is
   // otherwise not used.
   std::vector<std::pair<double, double> > timeCuts; 
   m_exposure = new Likelihood::LikeExposure(m_pars["binsz"], 
                                             m_pars["dcostheta"],
                                             timeCuts, gtis, zmax, zmin);
   std::string scFile = m_pars["scfile"];
   st_facilities::Util::file_ok(scFile);
   std::vector<std::string> scFiles;
   st_facilities::Util::resolve_fits_files(scFile, scFiles);
   std::vector<std::string>::const_iterator scIt = scFiles.begin();
   for ( ; scIt != scFiles.end(); scIt++) {
      st_facilities::Util::file_ok(*scIt);
      formatter.err() << "Working on file " << *scIt << std::endl;
      const tip::Table * scData = 
         tip::IFileSvc::instance().readTable(*scIt, m_pars["sctable"],
                                             filter.str());
      formatter.info(4) << "read " << scData->getNumRecords() 
                        << " rows" << std::endl;
      int chatter = m_pars["chatter"];
      bool print_output(true);
      if (chatter < 2) {
         print_output = false;
      }
      m_exposure->load(scData, print_output);
      delete scData;
   }

   if (m_exposure->numIntervals() == 0) {
      formatter.warn() << "WARNING: No intervals have been read in from "
                       << "the FT2 files that correspond to the FT1 data.\n"
                       << "All livetimes will be identically zero."
                       << std::endl;
   }
}
