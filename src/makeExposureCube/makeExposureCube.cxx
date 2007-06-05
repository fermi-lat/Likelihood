/** 
 * @file makeExposureCube.cxx
 * @brief Create an Exposure hypercube.
 * @author J. Chiang
 *
 *  $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/makeExposureCube/makeExposureCube.cxx,v 1.43 2007/04/30 04:00:40 jchiang Exp $
 */

#include <cstdlib>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "st_stream/StreamFormatter.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "st_facilities/Util.h"

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
   void getTimeBounds(const std::vector< std::pair<double, double> > & gtis,
                      const std::vector< std::pair<double, double> > & tranges,
                      double & tmin, double & tmax) {
      getTBounds(gtis, tmin, tmax);
      double t0(tmin), t1(tmax);
      getTBounds(tranges, t0, t1);
      if (t0 < tmin) {
         tmin = t0;
      }
      if (t1 > tmax) {
         tmax = t1;
      }
   }
}

/**
 * @class ExposureCube
 * @brief Class to encapsulate methods for creating an exposure
 * hypercube in (ra, dec, cos_theta) using the LikeExposure class.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/makeExposureCube/makeExposureCube.cxx,v 1.43 2007/04/30 04:00:40 jchiang Exp $
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
   m_exposure->write(output_file);
   std::auto_ptr<tip::Table> 
      table(tip::IFileSvc::instance().editTable(output_file, "Exposure"));
   m_roiCuts->writeDssTimeKeywords(table->getHeader());
   m_roiCuts->writeGtiExtension(output_file);
   
   tip::Header & header(table->getHeader());
   header["TSTART"].set(m_roiCuts->minTime());
   header["TSTOP"].set(m_roiCuts->maxTime());
   header.erase("TNULL1");
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

   std::vector<std::pair<double, double> > timeCuts;
   std::vector<std::pair<double, double> > gtis;
   m_roiCuts->getTimeCuts(timeCuts);
   m_roiCuts->getGtis(gtis);

   double tmin, tmax;
   ::getTimeBounds(gtis, timeCuts, tmin, tmax);
   static double maxIntervalSize(30);
   tmin -= 2.*maxIntervalSize;
   tmax += 2.*maxIntervalSize;
   std::ostringstream filter;
   filter << std::setprecision(20);
   filter << "(START >= " << tmin << ") && (STOP <= " << tmax << ")";
   formatter.info(4) << "applying filter: " << filter.str() << std::endl;

   m_exposure = new Likelihood::LikeExposure(m_pars["pixel_size"], 
                                             m_pars["cos_theta_step"],
                                             timeCuts, gtis);
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
      int chatter = m_pars["chatter"];
      bool print_output(true);
      if (chatter < 2) {
         print_output = false;
      }
      m_exposure->load(scData, print_output);
      delete scData;
   }
}
