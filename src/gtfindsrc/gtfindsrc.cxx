/**
 * @file gtfindsrc.cxx
 * @brief Use Nelder-Mead algorithm to fit for a point source location.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/gtfindsrc/gtfindsrc.cxx,v 1.8 2007/02/25 19:54:03 jchiang Exp $
 */

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <stdexcept>

#include "st_stream/StreamFormatter.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "st_facilities/Util.h"

#include "optimizers/Amoeba.h"
#include "optimizers/dArg.h"
#include "optimizers/Exception.h"
#include "optimizers/Optimizer.h"
#include "optimizers/OptimizerFactory.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/LogLike.h"
#include "Likelihood/Util.h"


using namespace Likelihood;

/**
 * @class findSrc
 * @brief Class for encapsulating methods for fitting for a source
 * position.
 *
 */

class findSrc : public st_app::StApp {
public:
   findSrc();
   virtual ~findSrc() throw() {
      try {
         delete m_testSrc;
         delete m_opt;
         delete m_logLike;
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
   PointSource * m_testSrc;
   std::vector<std::string> m_eventFiles;

   void promptForInputData();
   void readEventData();
   void readSrcModel();
   void identifyTarget();
   void selectOptimizer();
   double fitPosition(double step=0.1);
   double errEst(const std::vector< std::vector<double> > & testPoints) const;
   void setTestSource();
};

void findSrc::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

st_app::StAppFactory<findSrc> myAppFactory("gtfindsrc");

findSrc::findSrc() : st_app::StApp(), m_helper(0), 
                     m_pars(st_app::StApp::getParGroup("gtfindsrc")), 
                     m_logLike(0), m_opt(0), m_testSrc(0) {}

void findSrc::run() {
   promptForInputData();
   m_helper = new AppHelpers(&m_pars, "UNBINNED");

   std::string expcube_file = m_pars["exposure_cube_file"];
   if (expcube_file != "none" && expcube_file != "") {
      ExposureCube & expCube = 
         const_cast<ExposureCube &>(m_helper->observation().expCube());
      expCube.readExposureCube(expcube_file);
   }

   ResponseFunctions & respFuncs = 
      const_cast<ResponseFunctions &>(m_helper->observation().respFuncs());
   respFuncs.setEdispFlag(false);

   std::string exposureFile = m_pars["exposure_map_file"];
   std::string eventFile = m_pars["evfile"];
   std::string evtable = m_pars["evtable"];
   st_facilities::Util::file_ok(eventFile);
   st_facilities::Util::resolve_fits_files(eventFile, m_eventFiles);
   bool compareGtis(false);
   bool relyOnStreams(false);
   std::string respfunc = m_pars["rspfunc"];
   bool skipEventClassCuts(respfunc != "DSS");
   for (unsigned int i = 1; i < m_eventFiles.size(); i++) {
      AppHelpers::checkCuts(m_eventFiles[0], evtable, m_eventFiles[i],
                            evtable, compareGtis, relyOnStreams,
                            skipEventClassCuts);
   }
   compareGtis = true;
   if (exposureFile != "none" && exposureFile != "") {
      AppHelpers::checkCuts(m_eventFiles, evtable, exposureFile, "",
                            compareGtis, relyOnStreams,
                            skipEventClassCuts);
   }
   if (expcube_file != "none" && expcube_file != "") {
      AppHelpers::checkTimeCuts(m_eventFiles, evtable, expcube_file, 
                                "Exposure", compareGtis);
   }
   m_helper->setRoi();
   m_helper->readScData();
   m_helper->readExposureMap();

   m_logLike = new LogLike(m_helper->observation());
   readEventData();
   readSrcModel();
   identifyTarget();
   selectOptimizer();
   m_pars.Save();
   fitPosition();
   std::string target = m_pars["target_source"];
   bool clobber = m_pars["clobber"];
   if (clobber && m_logLike->getSource(target)) {
      std::string srcModelFile = m_pars["source_model_file"];
      m_logLike->writeXml(srcModelFile);
   }
}

void findSrc::promptForInputData() {
   m_pars.Prompt("evfile");
   m_pars.Prompt("scfile");
   m_pars.Prompt("outfile");
   m_pars.Prompt("rspfunc");
   m_pars.Prompt("exposure_cube_file");
   m_pars.Prompt("exposure_map_file");
   m_pars.Prompt("source_model_file");
   std::string outfile = m_pars["outfile"];
   bool clobber = m_pars["clobber"];
   if (!clobber && st_facilities::Util::fileExists(outfile)) {
      throw std::runtime_error("Output file, " + outfile 
                               + " exists and clobber is set to no.");
   }
}

void findSrc::readEventData() {
   std::vector<std::string> eventFiles;
   std::string evfile = m_pars["evfile"];
   Util::file_ok(evfile);
   Util::resolve_fits_files(evfile, eventFiles);
   std::vector<std::string>::const_iterator evIt = eventFiles.begin();
   for ( ; evIt != eventFiles.end(); evIt++) {
      Util::file_ok(*evIt);
      m_helper->observation().eventCont().getEvents(*evIt);
   }
}

void findSrc::readSrcModel() {
   std::string srcModelFile = m_pars["source_model_file"];
   if (srcModelFile != "" && srcModelFile != "none") {
      Util::file_ok(srcModelFile);
      st_stream::StreamFormatter formatter("findSrc", "readSrcModel", 2);
      formatter.info() << "Building source model from "
                       << srcModelFile << std::endl;
      m_logLike->readXml(srcModelFile, m_helper->funcFactory());
      m_logLike->computeEventResponses();
// evaluate log-likelihood using the input source model
      std::cout << "-log-likelihood of input source model: "
                << -m_logLike->value() << std::endl;
   }
}

void findSrc::identifyTarget() {
   m_pars.Prompt("target_source");
   std::string target = m_pars["target_source"];
   Source * candidateSrc(0);
   try {
      candidateSrc = m_logLike->deleteSource(target);
   } catch (optimizers::Exception & eObj) {
      st_stream::StreamFormatter formatter("findSrc", "identifyTarget", 2);
      std::string mode = m_pars["mode"];
      if (!myAppFactory.getGuiMode() && mode != "h") {
         formatter.info() << "Source '" << target 
                          << "' not found in source model." << std::endl;
         formatter.info() << "Enter coordinates for test source:" << std::endl;
      }
      setTestSource();
      return;
   }
   m_testSrc = dynamic_cast<PointSource *>(candidateSrc);
   if (m_testSrc == 0) {
      m_logLike->addSource(candidateSrc);
      throw std::runtime_error("Target source " + target + 
                               " is not a point source.");
   }
}

void findSrc::selectOptimizer() {
   m_pars.Prompt("optimizer");
   m_pars.Prompt("fit_tolerance");
   m_pars.Prompt("amoeba_tolerance");
   std::string optimizer = m_pars["optimizer"];
   m_opt = optimizers::OptimizerFactory::instance().create(optimizer, 
                                                           *m_logLike);
}

class LikeFunc : public optimizers::Functor {
public:
   LikeFunc(optimizers::Optimizer & opt, LogLike & logLike, 
            PointSource & testSrc, st_stream::StreamFormatter & formatter, 
            const std::string & coordSys="CEL",
            double tol=1e-5, bool findMin=true)
      : m_opt(opt), m_logLike(logLike), m_testSrc(testSrc),
        m_formatter(formatter), m_coordSys(coordSys), 
        m_tol(tol), m_findMin(findMin) {
      m_logLike.addSource(&m_testSrc);
   }
   virtual ~LikeFunc() {}
   virtual double operator()(std::vector<double> &coords) {
      m_logLike.deleteSource(m_testSrc.getName());
      if (m_coordSys == "CEL") {
         m_testSrc.setDir(coords[0], coords[1], true, false);
      } else {
         m_testSrc.setGalDir(coords[0], coords[1], true, false);
      }
      m_logLike.addSource(&m_testSrc);
      if (m_findMin) {
         m_opt.find_min(0, m_tol);
      }
      optimizers::dArg dummy(1.);
      double test_value = -m_logLike(dummy);
      m_formatter.info().precision(10);
      m_formatter.info() << coords[0] << "  "
                         << coords[1] << "  "
                         << test_value << std::endl;
      std::vector<double> point;
      point.push_back(coords[0]);
      point.push_back(coords[1]);
      point.push_back(test_value);
      m_testPoints.push_back(point);
      return test_value;
   }
   const std::vector< std::vector<double> > & testPoints() const {
      return m_testPoints;
   }
   // Sort from highest to lowest function value.
   void sortPoints() {
      std::stable_sort(m_testPoints.begin(), m_testPoints.end(), compare);
   }
private:
   optimizers::Optimizer & m_opt;
   LogLike & m_logLike;
   PointSource & m_testSrc;
   st_stream::StreamFormatter & m_formatter;
   std::string m_coordSys;
   double m_tol;
   bool m_findMin;

   std::vector< std::vector<double> > m_testPoints;

   // Want to sort from highest to lowest.
   static bool compare(const std::vector<double> & x,
                       const std::vector<double> & y) {
      return x.at(2) > y.at(2);
   }
};

double findSrc::fitPosition(double step) {
   std::vector<double> coords(2);
   bool use_lb = m_pars["use_lb"];
   std::string coordSys("CEL");
   coords[0] = m_testSrc->getDir().ra();
   coords[1] = m_testSrc->getDir().dec();
   if (use_lb) {
      coordSys = "GAL";
      coords[0] = m_testSrc->getDir().l();
      coords[1] = m_testSrc->getDir().b();
   }
   double tol = m_pars["fit_tolerance"];
   bool reopt = m_pars["reoptimize"];
   st_stream::StreamFormatter formatter("findSrc", "fitPosition", 2);
   LikeFunc func(*m_opt, *m_logLike, *m_testSrc, formatter, coordSys, 
                 tol, reopt);
   bool addstep;
   optimizers::Amoeba my_amoeba(func, coords, step, addstep=true);
   double pos_tol = m_pars["amoeba_tolerance"];
   double statValue = my_amoeba.findMin(coords, pos_tol);
   if (m_testSrc->getName() == "testSource") {
      m_logLike->deleteSource("testSource");
   }
   func.sortPoints();
   double pos_error(0);
   try {
      pos_error = errEst(func.testPoints());
   } catch (std::runtime_error & eObj) {
      formatter.info() << eObj.what() << std::endl;
   }
   std::string outfile = m_pars["outfile"];
   const std::vector< std::vector<double> > & testPoints(func.testPoints());
   if (outfile != "" && outfile != "none") {
      std::ofstream output(outfile.c_str());
      output << std::setprecision(10);
      for (size_t i = 0; i < testPoints.size(); i++) {
         output << testPoints.at(i).at(0) << "  "
                << testPoints.at(i).at(1) << "  "
                << testPoints.at(i).at(2) << "  "
                << pos_error << std::endl;
      }
      output.close();
   }
   formatter.info() << "Best fit position: "
                    << testPoints.back().at(0) << ",  "
                    << testPoints.back().at(1) << "\n"
                    << "Error circle radius: " << pos_error
                    << std::endl;
   return statValue;
}

double findSrc::
errEst(const std::vector< std::vector<double> > & testPoints) const {
   double Sx(0);
   double Sy(0);
   double Sxy(0);
   double Sxx(0);
   size_t npts(0);
   double deltaL(2);
   double Lmin(testPoints.back().at(2));
   astro::SkyDir minDir(testPoints.back().at(0), testPoints.back().at(1));
   for (size_t i = 0; i < testPoints.size(); i++) {
      const std::vector<double> & z(testPoints.at(i));
      if (z.at(2) < Lmin + deltaL) {
         double x(minDir.difference(astro::SkyDir(z.at(0), z.at(1))));
         x = x*x;
         Sx += x;
         Sy += z.at(2);
         Sxy += x*z.at(2);
         Sxx += x*x;
         npts++;
      }
   }
   double numerator(npts*Sxy - Sy*Sx);
   double denominator(npts*Sxx - Sx*Sx);
   if (denominator == 0 || numerator == 0 || numerator/denominator <= 0 ) {
      throw std::runtime_error("A reliable positional error estimate cannot "
                               "be made.\nPlease inspect the output file");
   }
   double AA(numerator/denominator);
   return 180./M_PI/std::sqrt(2.*AA);
}

void findSrc::setTestSource() {
   m_pars.Prompt("use_lb");
   m_pars.Prompt("lon0");
   m_pars.Prompt("lat0");
   double ra = m_pars["lon0"];
   double dec = m_pars["lat0"];
   if (m_pars["use_lb"]) {
      astro::SkyDir srcDir(m_pars["lon0"], m_pars["lat0"],
                           astro::SkyDir::GALACTIC);
      ra = srcDir.ra();
      dec = srcDir.dec();
   }
   m_testSrc = new PointSource(ra, dec, m_helper->observation());
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
   m_testSrc->setSpectrum(pl);
   m_testSrc->setName("testSource");
}
