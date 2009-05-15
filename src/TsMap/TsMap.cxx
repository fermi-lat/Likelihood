/**
 * @file TsMap.cxx
 * @brief Application for producing "test-statistic" maps.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/TsMap/TsMap.cxx,v 1.43 2009/03/23 23:29:12 jchiang Exp $
 */

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <sstream>

#include "facilities/commonUtilities.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Image.h"

#include "st_stream/StreamFormatter.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "st_facilities/Util.h"

#include "optimizers/dArg.h"
#include "optimizers/Optimizer.h"
#include "optimizers/OptimizerFactory.h"
#include "optimizers/Exception.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/LogLike.h"

using namespace Likelihood;

/**
 * @class TsMap
 * @brief Class for encapsulating methods for creating a test-statistic
 * map.
 *
 */

class TsMap : public st_app::StApp {
public:
   TsMap();
   virtual ~TsMap() throw() {
      try {
         delete m_opt;
         delete m_helper;
      } catch (std::exception & eObj) {
         std::cerr << eObj.what() << std::endl;
      } catch (...) {
      }
   }
   virtual void run();
   virtual void banner() const;
private:
   AppHelpers * m_helper;
   st_app::AppParGroup & m_pars;
   LogLike * m_logLike;
   optimizers::Optimizer * m_opt;
   st_stream::StreamFormatter * m_formatter;

   std::vector<std::string> m_eventFiles;
   std::vector<astro::SkyDir> m_dirs;
   std::vector<float> m_tsMap;
   std::string m_coordSys;
   std::vector<double> m_crpix;
   std::vector<double> m_crval;
   std::vector<double> m_cdelt;
   void readSrcModel();
   void readEventData();
   void selectOptimizer();
   void setGrid();
   void computeMap();
   void writeFitsFile();
   void setPointSourceSpectrum(PointSource &src);

   static std::string s_cvs_id;
};

st_app::StAppFactory<TsMap> myAppFactory("gttsmap");

TsMap::TsMap() 
   : st_app::StApp(), m_helper(0), 
     m_pars(st_app::StApp::getParGroup("gttsmap")),
     m_logLike(0), m_opt(0),
     m_formatter(new st_stream::StreamFormatter("gttsmap", "", 2)) {
   setVersion(s_cvs_id);
}

std::string TsMap::s_cvs_id("$Name:  $");

void TsMap::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

void TsMap::run() {
   m_pars.Prompt();
   m_pars.Save();

   m_helper = new AppHelpers(&m_pars, "UNBINNED");
   m_helper->checkOutputFile();
   std::string coordSys = m_pars["coordsys"];
   m_coordSys = coordSys;
   std::string expCubeFile = m_pars["expcube"];
   if (expCubeFile != "" && expCubeFile != "none") {
      m_helper->observation().expCube().readExposureCube(expCubeFile);
   }
   st_facilities::Util::file_ok(m_pars["evfile"]);
   st_facilities::Util::resolve_fits_files(m_pars["evfile"], m_eventFiles);
   std::string ev_table = m_pars["evtable"];
   bool compareGtis;
   bool relyOnStreams;
   std::string respfunc = m_pars["irfs"];
   bool skipEventClassCuts(respfunc != "DSS");
   for (size_t i(1); i < m_eventFiles.size(); i++) {
      AppHelpers::checkCuts(m_eventFiles[0], ev_table,
                            m_eventFiles[i], ev_table,
                            compareGtis=false, 
                            relyOnStreams=false, 
                            skipEventClassCuts);
   }
   m_helper->setRoi();
   m_helper->readScData();
   m_helper->readExposureMap();

   m_logLike = new LogLike(m_helper->observation());
   readEventData();
   readSrcModel();

   selectOptimizer();
   setGrid();

   computeMap();
   writeFitsFile();
}

void TsMap::readEventData() {
   std::vector<std::string>::const_iterator evFile(m_eventFiles.begin());
   for ( ; evFile != m_eventFiles.end(); ++evFile) {
      st_facilities::Util::file_ok(*evFile);
      m_logLike->getEvents(*evFile);
   }
}

void TsMap::readSrcModel() {
   std::string srcModelFile = m_pars["srcmdl"];
   if (srcModelFile != "" && srcModelFile != "none") {
      st_facilities::Util::file_ok(srcModelFile);
      m_logLike->readXml(srcModelFile, m_helper->funcFactory());
      m_logLike->computeEventResponses();
   }
}   

void TsMap::selectOptimizer() {
   std::string optimizer = m_pars["optimizer"];
   m_opt = optimizers::OptimizerFactory::instance().create(optimizer,
                                                           *m_logLike);
}

void TsMap::setGrid() {
   int nxpix = m_pars["nxpix"];
   int nypix = m_pars["nypix"];
   double xref = m_pars["xref"];
   double yref = m_pars["yref"];
   double binsize = m_pars["binsz"];
   std::string coordsys = m_pars["coordsys"];
   bool is_galactic(coordsys == "GAL");
   std::string proj_name = m_pars["proj"];

   double crpix[] = {nxpix/2. + 0.5, nypix/2. + 0.5};
   m_crpix = std::vector<double>(crpix, crpix+2);
   double crval[] = {xref, yref};
   m_crval = std::vector<double>(crval, crval+2);
   double cdelt[] = {binsize, binsize};
   m_cdelt = std::vector<double>(cdelt, cdelt+2);
   astro::SkyProj proj(proj_name, crpix, crval, cdelt, 0, is_galactic);

   m_dirs.clear();
   for (int j(0); j < nypix; j++) {
      for (int i(0); i < nxpix; i++) {
         m_dirs.push_back(astro::SkyDir(i+1, j+1, proj));
      }
   }
}

void TsMap::computeMap() {
   Likelihood::PointSource testSrc(m_dirs.at(0).ra(), m_dirs.at(0).dec(), 
                                   m_helper->observation());
   setPointSourceSpectrum(testSrc);
   testSrc.setName("testSource");

   int verbosity = m_pars["chatter"];
   verbosity -= 2;
   double tol = m_pars["ftol"];
   std::string tol_type = m_pars["toltype"];
   optimizers::TOLTYPE tolType(optimizers::ABSOLUTE);
   if (tol_type == "REL") {
      tolType = optimizers::RELATIVE;
   }
   double logLike0;
   try {
      m_opt->find_min_only(verbosity, tol, tolType);
      logLike0 = m_logLike->value();
   } catch (...) {
      logLike0 = 0;
   }
   bool computeExposure(true);

   int step(m_dirs.size()/20);
   if (step == 0) {
      step = 2;
   }

   for (size_t i(0); i < m_dirs.size(); i++) {
      if ((i % step) == 0) {
         m_formatter->warn() << ".";
      }
      testSrc.setDir(m_dirs.at(i).ra(), m_dirs.at(i).dec(),
                     computeExposure, false);

      m_logLike->addSource(&testSrc);
      try {
         m_opt->find_min_only(verbosity, tol, tolType);
         m_tsMap.push_back(2.*(m_logLike->value() - logLike0));
      } catch (optimizers::Exception & eObj) {
         m_formatter->err() << eObj.what() << std::endl;
         // Default null value.
         m_tsMap.push_back(0);
      }
      m_formatter->info(3) << m_dirs.at(i).ra() << "  "
                           << m_dirs.at(i).dec() << "  "
                           << m_tsMap.back() << std::endl;
      m_logLike->deleteSource(testSrc.getName());
   }
   m_formatter->warn() << "!" << std::endl;
}

void TsMap::setPointSourceSpectrum(PointSource &src) {
   optimizers::Function * pl = m_helper->funcFactory().create("PowerLaw");
   double parValues[] = {1., -2., 100.};
   std::vector<double> pars(parValues, parValues + 3);
   pl->setParamValues(pars);
   optimizers::Parameter indexParam = pl->getParam("Index");
   indexParam.setBounds(-3.5, -1.);
   pl->setParam(indexParam);
   optimizers::Parameter prefactorParam = pl->getParam("Prefactor");
   prefactorParam.setBounds(1e-10, 1e3);
   prefactorParam.setScale(1e-9);
   pl->setParam(prefactorParam);
   src.setSpectrum(pl);
}

void TsMap::writeFitsFile() {
   std::string outfile = m_pars["outfile"];
   if (st_facilities::Util::fileExists(outfile)) {
      std::remove(outfile.c_str());
   }
   std::string dataPath = 
      facilities::commonUtilities::getDataPath("Likelihood");
   std::string templateFile = 
      facilities::commonUtilities::joinPath(dataPath, "TsMapTemplate");

   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   fileSvc.createFile(outfile, templateFile);

   tip::Image * image(fileSvc.editImage(outfile, ""));

   std::vector<long> naxes;
   naxes.push_back(m_pars["nxpix"]);
   naxes.push_back(m_pars["nypix"]);
   image->setImageDimensions(naxes);

   tip::Header & header(image->getHeader());
   header["TELESCOP"].set("GLAST");
   header["INSTRUME"].set("LAT SIMULATION");

   header["DATE-OBS"].set("");
   header["DATE-END"].set("");


   header["CRVAL1"].set(m_crval.at(0));
   header["CRVAL2"].set(m_crval.at(1));

   header["CRPIX1"].set(m_crpix.at(0));
   header["CRPIX2"].set(m_crpix.at(1));

   header["CDELT1"].set(m_cdelt.at(0));
   header["CDELT2"].set(m_cdelt.at(1));

   header["CROTA2"].set(0);

   header["CREATOR"].set("gttsmap " + getVersion());

   std::string proj = m_pars["proj"];
   std::string coordsys = m_pars["coordsys"];
   if (coordsys == "GAL") {
      header["CTYPE1"].set("GLON-" + proj);
      header["CTYPE2"].set("GLAT-" + proj);
   } else {
      header["CTYPE1"].set("RA---" + proj);
      header["CTYPE2"].set("DEC--" + proj);
   }
   image->set(m_tsMap);
   double tstart(m_helper->observation().roiCuts().minTime());
   double tstop(m_helper->observation().roiCuts().maxTime());
   st_facilities::Util::writeDateKeywords(image, tstart, tstop, false);
   delete image;
}
