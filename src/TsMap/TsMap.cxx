/**
 * @file TsMap.cxx
 * @brief Application for producing "test-statistic" maps.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/TsMap/TsMap.cxx,v 1.53 2012/11/11 03:26:02 jchiang Exp $
 */

#include <cmath>
#include <cstdio>
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

#include "st_facilities/Environment.h"
#include "st_facilities/Util.h"

#include "optimizers/dArg.h"
#include "optimizers/Optimizer.h"
#include "optimizers/OptimizerFactory.h"
#include "optimizers/Exception.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/BinnedLikelihood.h"
#include "Likelihood/LogLike.h"
#include "Likelihood/SourceMap.h"

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
   std::string m_statistic;
   CountsMap * m_dataMap;
   std::vector<astro::SkyDir> m_dirs;
   std::vector<float> m_tsMap;
   std::string m_coordSys;
   std::vector<double> m_crpix;
   std::vector<double> m_crval;
   std::vector<double> m_cdelt;
   void promptForParameters();
   void readSrcModel();
   void readEventData(const std::vector<std::string> & evfiles);
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
   m_pars.setSwitch("statistic");
   m_pars.setCase("statistic", "BINNED", "cmap");
   m_pars.setCase("statistic", "BINNED", "bexpmap");
   m_pars.setCase("statistic", "BINNED", "psfcorr");
   m_pars.setCase("statistic", "UNBINNED", "evfile");
   m_pars.setCase("statistic", "UNBINNED", "evtable");
   m_pars.setCase("statistic", "UNBINNED", "scfile");
   m_pars.setCase("statistic", "UNBINNED", "sctable");
   m_pars.setCase("statistic", "UNBINNED", "expmap");
}

std::string TsMap::s_cvs_id("$Name: Likelihood-18-00-04 $");

void TsMap::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

void TsMap::promptForParameters() {
   m_pars.Prompt("statistic");
   std::string statistic = m_pars["statistic"];
   m_statistic = statistic;
   if (m_statistic == "BINNED") {
      m_pars.Prompt("cmap");
      m_pars.Prompt("bexpmap");
   } else {
      m_pars.Prompt("scfile");
      m_pars.Prompt("evfile");
      m_pars.Prompt("expmap");
   }
   m_pars.Prompt("expcube");
   m_pars.Prompt("srcmdl");
   m_pars.Prompt("irfs");
   m_pars.Prompt("optimizer");

   m_pars.Prompt("outfile");
   m_pars.Prompt("nxpix");
   m_pars.Prompt("nypix");
   m_pars.Prompt("binsz");
   m_pars.Prompt("coordsys");
   m_pars.Prompt("xref");
   m_pars.Prompt("yref");
   m_pars.Prompt("proj");

   m_pars.Save();
}

void TsMap::run() {
   promptForParameters();

   m_helper = new AppHelpers(&m_pars, m_statistic);

   m_helper->checkOutputFile();
   std::string expcube = m_pars["expcube"];
   std::string irfs = m_pars["irfs"];
   if (expcube != "" && expcube != "none") {
      m_helper->observation().expCube().readExposureCube(expcube);
   }
   if (m_statistic == "UNBINNED") {
      std::vector<std::string> evfiles;
      st_facilities::Util::resolve_fits_files(m_pars["evfile"], evfiles);
      std::string evtable = m_pars["evtable"];
      bool compareGtis;
      bool relyOnStreams;
      bool skipEventClassCuts(irfs != "DSS");
      for (size_t i(1); i < evfiles.size(); i++) {
         AppHelpers::checkCuts(evfiles[0], evtable,
                               evfiles[i], evtable,
                               compareGtis=false, 
                               relyOnStreams=false, 
                               skipEventClassCuts);
      }
      m_helper->setRoi();
      m_helper->readScData();
      m_helper->readExposureMap();
      m_logLike = new LogLike(m_helper->observation());
      readEventData(evfiles);
   } else { // Assume we are operating in binned mode.
      std::string cmap = m_pars["cmap"];
      m_helper->setRoi(cmap, "", false);
      if (!m_helper->observation().expCube().haveFile()) {
         throw std::runtime_error
            ("An exposure cube file is required for binned analysis. "
             "Please specify an exposure cube file.");
      }
      st_facilities::Util::file_ok(cmap);
      m_dataMap = new CountsMap(cmap);
      bool apply_psf_corrections = m_pars["psfcorr"];
      bool computePointSources(true);
      m_logLike = new BinnedLikelihood(*m_dataMap, 
                                       m_helper->observation(),
                                       cmap, 
                                       computePointSources,
                                       apply_psf_corrections);
      std::string bexpmap = m_pars["bexpmap"];
      AppHelpers::checkExposureMap(cmap, bexpmap);
      dynamic_cast<BinnedLikelihood *>(m_logLike)->setVerbose(false);
   }
   readSrcModel();
   selectOptimizer();
   setGrid();
   computeMap();
   writeFitsFile();
}

void TsMap::readEventData(const std::vector<std::string> & evfiles) {
   std::vector<std::string>::const_iterator evfile(evfiles.begin());
   for ( ; evfile != evfiles.end(); ++evfile) {
      st_facilities::Util::file_ok(*evfile);
      m_logLike->getEvents(*evfile);
   }
}

void TsMap::readSrcModel() {
   std::string srcModelFile = m_pars["srcmdl"];
   if (srcModelFile != "" && srcModelFile != "none") {
      st_facilities::Util::file_ok(srcModelFile);
      bool requireExposure = (m_statistic != "BINNED");
//      bool loadMaps = (m_statistic != "BINNED");
      bool loadMaps;
      bool addPointSources;
      m_logLike->readXml(srcModelFile, m_helper->funcFactory(),
                         requireExposure, addPointSources=true,
                         loadMaps=false);
      if (m_statistic == "UNBINNED") {
         m_logLike->computeEventResponses();
      }
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
   double cdelt[] = {-binsize, binsize};
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
   Likelihood::PointSource * testSrc(0);
   if (m_statistic == "UNBINNED") {
      testSrc = new Likelihood::PointSource(m_dirs.at(0).ra(), 
                                            m_dirs.at(0).dec(), 
                                            m_helper->observation());
   } else {
      testSrc = new Likelihood::PointSource();
   }
   setPointSourceSpectrum(*testSrc);
   testSrc->setName("testSource");

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
   int step(m_dirs.size()/20);
   if (step == 0) {
      step = 2;
   }
   bool computeExposure;
   for (size_t i(0); i < m_dirs.size(); i++) {
      if ((i % step) == 0) {
         m_formatter->warn() << ".";
      }
      testSrc->setDir(m_dirs.at(i).ra(), m_dirs.at(i).dec(),
                      computeExposure=(m_statistic=="UNBINNED"), false);

      m_logLike->addSource(testSrc);
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
      m_logLike->deleteSource(testSrc->getName());
      if (m_statistic == "BINNED") {
         dynamic_cast<BinnedLikelihood *>(m_logLike)
            ->eraseSourceMap(testSrc->getName());
      }
   }
   m_formatter->warn() << "!" << std::endl;
   delete testSrc;
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
   std::string dataPath(st_facilities::Environment::dataPath("Likelihood"));
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
