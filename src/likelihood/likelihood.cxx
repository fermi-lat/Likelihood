/**
 * @file likelihood.cxx
 * @brief Prototype standalone application for the Likelihood tool.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/likelihood/likelihood.cxx,v 1.160 2012/04/14 20:59:08 jchiang Exp $
 */

#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>

#include "st_stream/StreamFormatter.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "st_facilities/Util.h"

#include "st_graph/Engine.h"
#include "st_graph/IFrame.h"
#include "st_graph/Placer.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "optimizers/Optimizer.h"
#include "optimizers/OptimizerFactory.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/BinnedLikelihood.h"
#include "Likelihood/BinnedLikelihood2.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/CountsSpectra.h"
#include "Likelihood/ExposureCube.h"
#include "Likelihood/LogLike.h"
#include "Likelihood/OptEM.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/Source.h"
#include "Likelihood/SourceMap.h"

#include "EasyPlot.h"
#include "MathUtil.h"

namespace {
   double ptsrcSeparation(Likelihood::Source * src1,
                          Likelihood::Source * src2) {
      Likelihood::PointSource * 
         ptsrc1(dynamic_cast<Likelihood::PointSource *>(src1));
      Likelihood::PointSource * 
         ptsrc2(dynamic_cast<Likelihood::PointSource *>(src2));
      if (ptsrc1 == 0  || ptsrc2 == 0) {
         throw std::runtime_error("likelihood::ptsrcSeparation: Attempt "
                                  "to compute separation between two Source "
                                  "objects that are not both PointSources.");
      }
      return ptsrc1->getDir().difference(ptsrc2->getDir())*180./M_PI;
   }
   void strip_at_sign(std::string & input) {
      if (input.find_first_of("@") == 0) {
         std::string output = "";
         std::string::iterator it = input.begin() + 1;
         for ( ; it != input.end(); ++it) {
            output += *it;
         }
         input = output;
      }
   }
}

using namespace Likelihood;

/**
 * @class likelihood
 * @brief A class encapsulating the methods for performing a
 * Likelihood analysis in ballistic fashion.
 *
 * @author J. Chiang
 */

class likelihood : public st_app::StApp {
public:
   likelihood();
   virtual ~likelihood() throw() {
      try {
         delete m_logLike;
         delete m_opt;
         delete m_dataMap;
         delete m_formatter;
      } catch (std::exception & eObj) {
         std::cout << eObj.what() << std::endl;
      } catch (...) {
      }
   }
   virtual void run();
   virtual void banner() const;
   double cputime() const {
      return static_cast<double>(std::clock() - m_cpuStart)/1e6;
   }

private:
   AppHelpers * m_helper;
   st_app::AppParGroup & m_pars;

   LogLike * m_logLike;
   optimizers::Optimizer * m_opt;
   std::vector<std::string> m_eventFiles;
   CountsMap * m_dataMap;

   st_stream::StreamFormatter * m_formatter;

   std::string m_statistic;

   std::clock_t m_cpuStart;

   double m_tol;

   optimizers::TOLTYPE m_tolType;

   std::vector< std::vector<double> > m_covarianceMatrix;

   void promptForParameters();
   void createStatistic();
   void readEventData();
   void readSourceModel();
   void selectOptimizer(std::string optimizer="");
   void writeSourceXml();
   void writeCountsSpectra();
   void plotCountsSpectra();
   void printFitResults(const std::vector<double> &errors);
   void printFitQuality() const;
   bool prompt(const std::string &query);
   void setErrors(const std::vector<double> & errors);

   void computeTsValues(const std::vector<std::string> & srcNames,
                        std::map<std::string, double> & TsValues, 
                        std::map<std::string, double> & RoiDist);
   Source * m_tsSrc;
   double m_maxdist;
   void renormModel();
   void npredValues(double &, double &) const;
   optimizers::Parameter & normPar(Source *) const;
   bool isDiffuseOrNearby(Source *) const;
   double observedCounts();
   double sourceFlux(const std::string & srcName, double & fluxError);

   void getElims(double & emin, double & emax) const;

   static std::string s_cvs_id;
};

st_app::StAppFactory<likelihood> myAppFactory("gtlike");

std::string likelihood::s_cvs_id("$Name:  $");

void likelihood::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

likelihood::likelihood() 
   : st_app::StApp(), m_helper(0), 
     m_pars(st_app::StApp::getParGroup("gtlike")),
     m_logLike(0), m_opt(0), m_dataMap(0), 
     m_formatter(new st_stream::StreamFormatter("gtlike", "", 2)),
     m_cpuStart(std::clock()), m_tolType(optimizers::RELATIVE),
     m_tsSrc(0), m_maxdist(20.) {
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

void likelihood::run() {
   promptForParameters();
   std::string statistic = m_pars["statistic"];
   std::string tolType = m_pars["toltype"];
   if (tolType == "REL") {
      m_tolType = optimizers::RELATIVE;
   } else {
      m_tolType = optimizers::ABSOLUTE;
   }

   m_helper = new AppHelpers(&m_pars, statistic);
   std::string expcube_file = m_pars["expcube"];
   // if (expcube_file != "none" && expcube_file != "") {
   //    ExposureCube & expCube = 
   //       const_cast<ExposureCube &>(m_helper->observation().expCube());
   //    expCube.readExposureCube(expcube_file);
   // }
//   bool useEdisp = m_pars["edisp"];
   bool useEdisp = false;
   ResponseFunctions & respFuncs = 
      const_cast<ResponseFunctions &>(m_helper->observation().respFuncs());
   respFuncs.setEdispFlag(useEdisp);
   if (m_statistic == "BINNED") {
      m_helper->setRoi(m_pars["cmap"], "", false);
   } else {
      std::string exposureFile = m_pars["expmap"];
      std::string eventFile = m_pars["evfile"];
      std::string evtable = m_pars["evtable"];
      ::strip_at_sign(eventFile);
      st_facilities::Util::file_ok(eventFile);
      st_facilities::Util::resolve_fits_files(eventFile, m_eventFiles);
      bool compareGtis(false);
      bool relyOnStreams(false);
      std::string respfunc = m_pars["irfs"];
      bool skipEventClassCuts(respfunc != "DSS");
      for (size_t i = 1; i < m_eventFiles.size(); i++) {
         AppHelpers::checkCuts(m_eventFiles[0], evtable, m_eventFiles[i],
                               evtable, compareGtis, relyOnStreams,
                               skipEventClassCuts);
      }
      compareGtis = true;
      if (exposureFile != "none" && exposureFile != "") {
         AppHelpers::checkExpMapCuts(m_eventFiles, exposureFile, evtable, "");
      }
      if (expcube_file != "none" && expcube_file != "") {
         AppHelpers::checkTimeCuts(m_eventFiles, evtable, expcube_file, 
                                   "Exposure", compareGtis);
      }
      m_helper->setRoi();
      m_helper->readScData();
      m_helper->readExposureMap();
   }
   createStatistic();

// Set the verbosity level and convergence tolerance.
   long verbose = m_pars["chatter"];
// ST proclaimed nominal verbosity level is 2, but optimizers expect 1, so
// we subtract 1.
   if (verbose > 1) verbose--;
   double tol = m_pars["ftol"];
   m_tol = tol;
   std::vector<double> errors;

// The fit loop.  If indicated, query the user at the end of each
// iteration whether the fit is to be performed again.  This allows
// the user to adjust the source model xml file by hand between
// iterations.
   bool queryLoop = m_pars["refit"];
   do {
      readSourceModel();
// Do the fit.
/// @todo Allow the optimizer to be re-selected here by the user.    
      selectOptimizer();
      try {
         m_opt->find_min(verbose, m_tol, m_tolType);
         try {
            errors = m_opt->getUncertainty();
            setErrors(errors);
         } catch (optimizers::Exception & eObj) {
            m_formatter->err() << "Exception encountered while "
                               << "estimating errors:\n"
                               << eObj.what() << std::endl;
//             throw;
         }
      } catch (optimizers::Exception & eObj) {
         m_formatter->err() << "Exception encountered while minimizing "
                            << "objective function:\n"
                            << eObj.what() << std::endl;
//          throw;
      }
      try {
         m_covarianceMatrix = m_opt->covarianceMatrix();
      } catch (std::runtime_error &) {
         m_covarianceMatrix.clear();
      }

      printFitResults(errors);
      writeSourceXml();
      if (m_pars["plot"]) {
         plotCountsSpectra();
      }
   } while (queryLoop && prompt("Refit? [y] "));
   if (m_pars["save"]) {
      writeCountsSpectra();
   }
   m_formatter->info() << "Elapsed CPU time: " << cputime() << std::endl;
   delete m_helper;
}

void likelihood::setErrors(const std::vector<double> & errors) {
   std::vector<optimizers::Parameter> params;
   m_logLike->getFreeParams(params);
   if (errors.size() != params.size()) {
      throw std::runtime_error("number of error estimates does not "
                               "match the number of free parameters");
   }
   for (unsigned int i = 0; i < errors.size(); i++) {
      params.at(i).setError(errors.at(i));
   }
   m_logLike->setFreeParams(params);
}

void likelihood::promptForParameters() {
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
   m_pars.Prompt("sfile");
   AppHelpers::checkOutputFile(m_pars["clobber"], 
                               m_pars["sfile"]);
//    m_pars.Prompt("fluxmdl");
//    AppHelpers::checkOutputFile(m_pars["clobber"], 
//                                m_pars["fluxmdl"]);
   m_pars.Prompt("irfs");
//   m_pars.Prompt("edisp");
   m_pars.Prompt("optimizer");
   m_pars.Prompt("save");
   m_pars.Prompt("refit");
   m_pars.Prompt("toltype");

   m_pars.Save();
}

void likelihood::createStatistic() {
   if (m_statistic == "BINNED") {
      if (!m_helper->observation().expCube().haveFile()) {
         throw std::runtime_error
            ("An exposure cube file is required for binned analysis. "
             "Please specify an exposure cube file.");
      }
      std::string countsMapFile = m_pars["cmap"];
      st_facilities::Util::file_ok(countsMapFile);
      m_dataMap = new CountsMap(countsMapFile);
      bool apply_psf_corrections(false);
      bool computePointSources(true);
      try {
         apply_psf_corrections = m_pars["psfcorr"];
      } catch (...) {
         // assume parameter does not exist, so use default value.
      }
      if (::getenv("USE_BINNED_LIKELIHOOD2")) {
         m_logLike = new BinnedLikelihood2(*m_dataMap, m_helper->observation(),
                                           countsMapFile, 
                                           computePointSources,
                                           apply_psf_corrections);
      } else {
         m_logLike = new BinnedLikelihood(*m_dataMap, m_helper->observation(),
                                          countsMapFile, 
                                          computePointSources,
                                          apply_psf_corrections);
      }
      std::string binnedMap = m_pars["bexpmap"];
      AppHelpers::checkExposureMap(m_pars["cmap"], m_pars["bexpmap"]);
      return;
   } else if (m_statistic == "UNBINNED") {
      m_logLike = new LogLike(m_helper->observation());
   }
   readEventData();
}

void likelihood::readEventData() {
   std::vector<std::string>::const_iterator evIt = m_eventFiles.begin();
   for ( ; evIt != m_eventFiles.end(); evIt++) {
      st_facilities::Util::file_ok(*evIt);
      if (m_statistic == "BINNED")
	m_helper->observation().eventCont().getEvents(*evIt); // Why?
      else if (m_statistic == "UNBINNED")
	m_logLike->getEvents(*evIt);
   }
}

void likelihood::readSourceModel() {
   std::string sourceModel = m_pars["srcmdl"];
   bool requireExposure(true);
   bool addPointSources;
   bool loadMaps;
   if (m_statistic == "BINNED") {
      requireExposure = false;
   }
   if (m_logLike->getNumSrcs() == 0) {
// Read in the Source model for the first time.
      st_facilities::Util::file_ok(sourceModel);
      m_logLike->readXml(sourceModel, m_helper->funcFactory(), 
                         requireExposure, addPointSources=true,
                         loadMaps=false);
//                         loadMaps=(m_statistic!="BINNED"));
      if (m_statistic != "BINNED") {
         m_logLike->computeEventResponses();
      }
   } else {
// Re-read the Source model from the xml file, allowing only for 
// Parameter adjustments.
      st_facilities::Util::file_ok(sourceModel);
      m_logLike->reReadXml(sourceModel);
   }
}

void likelihood::selectOptimizer(std::string optimizer) {
   delete m_opt;
   m_opt = 0;
   if (optimizer == "") {
// Type conversions for hoops do not work properly, so we are
// forced to use an intermediate variable.
      std::string opt = m_pars["optimizer"];
      optimizer = opt;
   }
   m_opt = optimizers::OptimizerFactory::instance().create(optimizer,
                                                           *m_logLike);
}

void likelihood::writeSourceXml() {
   std::string xmlFile = m_pars["sfile"];
   std::string funcFileName("");
   if (xmlFile != "none") {
      m_formatter->info() << "Writing fitted model to " 
                          << xmlFile << std::endl;
      m_logLike->writeXml(xmlFile, funcFileName);
   }
}

class EventData {
public:
   EventData(const std::vector<Event> & events) : m_events(events) {}
   unsigned long nobs(double emin, double emax) const {
      unsigned long nevents(0);
      std::vector<Event>::const_iterator it = m_events.begin();
      for ( ; it != m_events.end(); ++it) {
         if (emin <= it->getEnergy() && it->getEnergy() <= emax) {
            nevents++;
         }
      }
      return nevents;
   }
private:
   const std::vector<Event> & m_events;
};

void likelihood::writeCountsSpectra() {
   CountsSpectra counts(*m_logLike);

   if (m_statistic == "UNBINNED") {
      const RoiCuts & roiCuts = m_helper->observation().roiCuts();
      double emin(roiCuts.getEnergyCuts().first);
      double emax(roiCuts.getEnergyCuts().second);
      counts.setEbounds(emin, emax, 21);
   }

   std::string outfile = m_pars["specfile"];
   counts.writeTable(outfile);
}

void likelihood::plotCountsSpectra() {
   CountsSpectra counts(*m_logLike);

   if (m_statistic == "UNBINNED") {
      const RoiCuts & roiCuts = m_helper->observation().roiCuts();
      double emin(roiCuts.getEnergyCuts().first);
      double emax(roiCuts.getEnergyCuts().second);
      counts.setEbounds(emin, emax, 21);
   }
   
   std::vector<double> nobs;
   counts.getObsCounts(nobs);
   std::vector<double> nobs_err;

   std::vector<double> ebounds(counts.ebounds());
   std::vector<double> evals;
   std::vector<double> ewidth;
   for (size_t k = 0; k < ebounds.size()-1; k++) {
      evals.push_back(std::sqrt(ebounds.at(k)*ebounds.at(k+1)));
      ewidth.push_back(ebounds.at(k+1)-ebounds.at(k));
      nobs_err.push_back(std::sqrt(nobs.at(k)));
   }

// Normalize counts by bin-width to get counts/energy.
   for (size_t k = 0; k < ewidth.size(); k++) {
      nobs[k] /= ewidth[k];
      nobs_err[k] /= ewidth[k];
   }
   
// Set up a finer spectrum, in order to plot the predicted counts
// with higher resolution.
   CountsSpectra fine_counts(*m_logLike);
   std::vector<double>::size_type num_fine_points = ebounds.size() * 1;
   if (m_statistic == "UNBINNED") {
      fine_counts.setEbounds(ebounds.front(), ebounds.back(), num_fine_points);
   }

   std::vector<double> fine_ebounds(fine_counts.ebounds());
   std::vector<double> fine_evals;
   std::vector<double> fine_ewidth;
   for (size_t k = 0; k < fine_ebounds.size()-1; k++) {
      fine_evals.push_back(std::sqrt(fine_ebounds.at(k)*fine_ebounds.at(k+1)));
      fine_ewidth.push_back(fine_ebounds.at(k+1)-fine_ebounds.at(k));
   }

   std::vector<std::string> sourceNames;
   m_logLike->getSrcNames(sourceNames);
   std::vector<std::vector<double> > npred(sourceNames.size());
   std::vector<std::vector<double> > fine_npred(sourceNames.size());
   for (size_t i = 0; i < sourceNames.size(); i++) {
      counts.getSrcCounts(sourceNames.at(i), npred.at(i));
      fine_counts.getSrcCounts(sourceNames.at(i), fine_npred.at(i));
   }

   std::vector<double> zero(fine_evals.size(), 1.);

   try {
// Create the main frame to hold all plot frames.
      st_graph::Engine & engine(st_graph::Engine::instance());
      std::auto_ptr<st_graph::IFrame> 
         mainFrame(engine.createMainFrame(0, 600, 600));
// Create the plot.
      EasyPlot plot(mainFrame.get(), "", true, true, "Energy (MeV)", 
                    "counts/MeV", 600, 400);
      plot.scatter(evals, nobs, nobs_err);
      std::vector<double> npred_tot(npred.at(0).size(), 0);
      std::vector<double> fine_npred_tot(fine_npred.at(0).size(), 0);
      int color = st_graph::Color::eBlack;
      for (unsigned int i = 0; i < npred.size(); i++) {
         for (size_t k = 0; k < ewidth.size(); k++) {
            npred[i][k] /= ewidth[k];
         }
         for (size_t k = 0; k < fine_ewidth.size(); k++) {
            fine_npred[i][k] /= fine_ewidth[k];
         }
         color = st_graph::Color::nextColor(color);
         plot.linePlot(fine_evals, fine_npred[i], color);
         for (size_t k = 0; k < npred_tot.size(); k++) {
            npred_tot.at(k) += npred.at(i).at(k);
         }
         for (size_t k = 0; k < fine_npred_tot.size(); k++) {
            fine_npred_tot.at(k) += fine_npred.at(i).at(k);
         }
      }
      plot.linePlot(fine_evals, fine_npred_tot);
      
      EasyPlot residuals_plot(mainFrame.get(), "", true, false, 
                              "Energy (MeV)", "(counts - model) / model",
                              600, 200);
      std::vector<double> zero(evals.size(), 0.);

      st_graph::TopEdge top_edge(residuals_plot.getPlotFrame());
      top_edge.below(st_graph::BottomEdge(plot.getPlotFrame()));

      std::vector<double> residuals(npred_tot.size(), 0.);
      std::vector<double> residuals_err(npred_tot.size(), 0.);
      for (unsigned int i = 0; i < npred_tot.size(); ++i) {
         residuals.at(i) = (nobs.at(i) - npred_tot.at(i))/npred_tot.at(i);
         residuals_err.at(i) = nobs_err.at(i)/npred_tot.at(i);
      }

      residuals_plot.scatter(evals, residuals, residuals_err);
      residuals_plot.linePlot(evals, zero, st_graph::Color::eBlack, "dashed");
      
      EasyPlot::run();
   } catch (std::exception &eObj) {
      std::string message = "RootEngine could not create";
      if (!st_facilities::Util::expectedException(eObj, message)) {
         throw;
      }
   }
}

void likelihood::getElims(double & emin, double & emax) const {
   if (m_statistic == "UNBINNED") {
      const RoiCuts & roiCuts = m_helper->observation().roiCuts();
      emin = roiCuts.getEnergyCuts().first;
      emax = roiCuts.getEnergyCuts().second;
   } else {
      CountsSpectra counts(*m_logLike);
      emin = counts.ebounds().front();
      emax = counts.ebounds().back();
   }
}

double likelihood::sourceFlux(const std::string & srcName, double & fluxError) {
   double emin, emax;
   getElims(emin, emax);

   Source * my_source(m_logLike->getSource(srcName));

   if (my_source->spectrum().getNumFreeParams() == 0 ||
       m_covarianceMatrix.empty()) {
      fluxError = 0;
      return my_source->flux(emin, emax);
   }

   std::vector<std::string> sourceNames;
   m_logLike->getSrcNames(sourceNames);

   std::map<std::string, int> par_index_map;
   int indx(0);
   for (size_t i(0); i < sourceNames.size(); i++) {
      const std::string & src(sourceNames.at(i));
      std::vector<std::string> parNames;
      m_logLike->getSource(src)->spectrum().getFreeParamNames(parNames);
      for (size_t j(0); j < parNames.size(); j++) {
         par_index_map[src + "::" + parNames.at(j)] = indx;
         indx++;
      }
   }
   if (m_covarianceMatrix.size() != par_index_map.size()) {
      throw std::runtime_error("Covariance matrix size does not match the "
                               "number of free parameters.");
   }

   std::vector< std::vector<double> > my_covar;
   std::vector<std::string> srcpars;
   my_source->spectrum().getFreeParamNames(srcpars);
   std::vector<std::string> pars;
   for (size_t i(0); i < srcpars.size(); i++) {
      pars.push_back(srcName + "::" + srcpars.at(i));
   }
   for (size_t i(0); i < pars.size(); i++) {
      int ix = par_index_map[pars.at(i)];
      std::vector<double> row;
      for (size_t j(0); j < pars.size(); j++) {
         row.push_back(m_covarianceMatrix.at(ix).at(par_index_map[pars.at(j)]));
      }
      my_covar.push_back(row);
   }
   std::vector<double> partials;
   size_t npts(1000);
   for (size_t i(0); i < srcpars.size(); i++) {
      partials.push_back(my_source->fluxDeriv(srcpars.at(i), emin, emax, npts));
   }
   std::vector<double> covar_dot_partials;
   for (size_t i(0); i < my_covar.size(); i++) {
      double value(0);
      for (size_t j(0); j < partials.size(); j++) {
         value += my_covar.at(i).at(j)*partials.at(j);
      }
      covar_dot_partials.push_back(value);
   }
   fluxError = 0;
   for (size_t i(0); i < partials.size(); i++) {
      fluxError += partials.at(i)*covar_dot_partials.at(i);
   }
   fluxError = std::sqrt(fluxError);
   return my_source->flux(emin, emax);
}

void likelihood::printFitResults(const std::vector<double> &errors) {
   std::vector<std::string> srcNames;
   m_logLike->getSrcNames(srcNames);
   
   std::map<std::string, double> TsValues;
   std::map<std::string, double> RoiDist;

   computeTsValues(srcNames, TsValues, RoiDist);

   std::vector<optimizers::Parameter> parameters;
   std::vector<double>::const_iterator errIt = errors.begin();

   std::string outfile = m_pars["results"];

   std::ofstream resultsFile(outfile.c_str());
   bool write_output_files = m_pars["save"];
   if (!write_output_files) {
      resultsFile.clear(std::ios::failbit);
   }

   resultsFile << "{";

   double emin, emax;
   getElims(emin, emax);

   m_formatter->info() << "\nPhoton fluxes are computed for the energy range " 
                       << emin << " to " << emax << " MeV" << std::endl;

   double totalNpred(0);
   for (unsigned int i = 0; i < srcNames.size(); i++) {
      Source * src = m_logLike->getSource(srcNames[i]);
      Source::FuncMap srcFuncs = src->getSrcFuncs();
      srcFuncs["Spectrum"]->getParams(parameters);
      m_formatter->info() << "\n" << srcNames[i] << ":\n";
      resultsFile << "'" << srcNames[i] << "': {";
      for (unsigned int j = 0; j < parameters.size(); j++) {
         m_formatter->info() << parameters[j].getName() << ": "
                             << parameters[j].getValue();
         resultsFile << "'" << parameters[j].getName() << "': "
                     << "'" << parameters[j].getValue();
         if (parameters[j].isFree() && errIt != errors.end()) {
            m_formatter->info() << " +/- " << *errIt;
            resultsFile << " +/- " << *errIt << "',\n";
            errIt++;
         } else {
            resultsFile << "',\n";
         }
         m_formatter->info() << std::endl;
      }
      if (m_statistic != "BINNED") {
         double npred(src->Npred());
         m_formatter->info() << "Npred: "
                             << npred << std::endl;
         resultsFile << "'Npred': '" << npred << "',\n";
         totalNpred += npred;
      }
      if (RoiDist.count(srcNames[i])) {
         m_formatter->info() << "ROI distance: "
                             << RoiDist[srcNames[i]] << std::endl;
         resultsFile << "'ROI distance': '" << RoiDist[srcNames[i]] << "',\n";
      }
      if (TsValues.count(srcNames[i])) {
         m_formatter->info() << "TS value: "
                             << TsValues[srcNames[i]] << std::endl;
         resultsFile << "'TS value': '" << TsValues[srcNames[i]] << "',\n";
      }
      double fluxError;
      double flux = sourceFlux(srcNames[i], fluxError);
      m_formatter->info() << "Flux: " << flux;
      resultsFile << "'Flux': '" << flux;
      if (fluxError > 0) {
         m_formatter->info() << " +/- " << fluxError;
         resultsFile << " +/- " << fluxError;
      }
      m_formatter->info() << " photons/cm^2/s" << std::endl;
      resultsFile << "',\n";
      resultsFile << "},\n";
   }
   resultsFile << "}" << std::endl;

   bool check_fit = m_pars["check_fit"];
   if (check_fit) {
// Check quality of the fit and issue warning if it appears not to be good.
      m_formatter->info() << "\n";
      printFitQuality();
   }

   m_formatter->info() << "\nTotal number of observed counts: "
                       << observedCounts() << "\n"
                       << "Total number of model events: ";
   if (m_statistic == "BINNED") {
      if (::getenv("USE_BINNED_LIKELIHOOD2")) {
         m_formatter->info() 
            << dynamic_cast<BinnedLikelihood2 *>(m_logLike)->npred();
      } else {
         m_formatter->info() 
            << dynamic_cast<BinnedLikelihood *>(m_logLike)->npred();
      }
   } else {
      m_formatter->info() << totalNpred;
   }
   m_formatter->info() << std::endl;
   resultsFile.close();

   if (!write_output_files) {
      std::remove(outfile.c_str());
   }

   m_formatter->info().precision(10);
   m_formatter->info() << "\n-log(Likelihood): "
//                       << std::setprecision(10)
                       << -m_logLike->value()
                       << "\n" << std::endl;
   delete m_opt;
   m_opt = 0;
}

void likelihood::printFitQuality() const {
   CountsSpectra countsSpec(*m_logLike);

   if (m_statistic == "UNBINNED") {
      const RoiCuts & roiCuts = m_helper->observation().roiCuts();
      double emin(roiCuts.getEnergyCuts().first);
      double emax(roiCuts.getEnergyCuts().second);
      countsSpec.setEbounds(emin, emax, 21);
   }

   std::vector<double> counts;
   std::vector<double> srcCounts;
   const std::vector<double> & ebounds(countsSpec.ebounds());
   countsSpec.getObsCounts(counts);
   countsSpec.getTotalSrcCounts(srcCounts);

   int lastEnergy = srcCounts.size();
   std::vector<double> & mean(srcCounts);
   std::vector<double> significance_val(lastEnergy, 0.);
   std::vector<double> threshold(lastEnergy, 0.);

   try {
      for (int k = 0; k < lastEnergy; k++) {
         if (counts[k]>0) {//counts>0
            threshold[k]=0.05;
            if (mean[k]>0) { //counts> 0 mean>0
               significance_val[k]=MathUtil::poissonSig(counts[k],mean[k]);
            } else {
               significance_val[k]=MathUtil::poissonSig(counts[k],mean[k]);
            }
         } else {//counts=0    
            if (mean[k]>0){
               significance_val[k]=1./mean[k]; 
               threshold[k]=0.5;//this number may change
            } else {//counts=mean=0
               significance_val[k]=1;//always >threshold
               threshold[k]=0.05;
            }
         }
      }
   
      bool aboveThreshold_sig = false;
      for (int k = 0; k < lastEnergy; k++) {
         if (significance_val[k] < threshold[k]) {
            if (!aboveThreshold_sig) {
               m_formatter->warn() << "WARNING: Fit may be bad in range [" 
                                   << ebounds.at(k) << ", ";
            }
            aboveThreshold_sig = true;
         } else if (aboveThreshold_sig) {
            aboveThreshold_sig = false;
            m_formatter->warn() << ebounds.at(k) << "] (MeV)" << std::endl;
         }
      }

     // Print the rest of the warning if the last point was above the threshold.
     if (aboveThreshold_sig)
        m_formatter->warn() << ebounds.at(lastEnergy) << "] (MeV)" << std::endl;

   } catch (const std::exception & x) {
      m_formatter->warn() << "Failed to compute Poisson significance: " << x.what() << std::endl;
   }
}

void likelihood::computeTsValues(const std::vector<std::string> & srcNames,
                                 std::map<std::string, double> & TsValues,
                                 std::map<std::string, double> & RoiDist) {
// Save current set of parameters.
   std::vector<double> fitParams;
   m_logLike->getFreeParamValues(fitParams);
   double logLike_value = m_logLike->value();

   int verbose(0);
   m_formatter->info() << "Computing TS values for each source ("
                       << srcNames.size() << " total)\n";
   astro::SkyDir roiCenter 
      = m_helper->observation().roiCuts().extractionRegion().center();
   for (unsigned int i = 0; i < srcNames.size(); i++) {
      m_formatter->warn() << ".";
      Source * src = m_logLike->getSource(srcNames[i]);
      if (src->getType() == "Point" &&
          src->spectrum().getNumFreeParams() > 0) {
         m_tsSrc = m_logLike->deleteSource(srcNames[i]);
         if (m_statistic != "BINNED") {
            RoiDist[srcNames[i]] = 
               dynamic_cast<PointSource *>(m_tsSrc)->getDir().
               difference(roiCenter)*180./M_PI;
         }
         if (m_logLike->getNumFreeParams() > 0) {
            selectOptimizer();
// Save value for null hypothesis before modifying parameters.
            double null_value(m_logLike->value());
            if (m_pars["tsmin"]) {
               try {
                  m_opt->find_min_only(verbose, m_tol, m_tolType);
               } catch (std::exception & eObj) {
                  m_formatter->err() << eObj.what() << std::endl;
               }
            } else {
               renormModel();
               m_logLike->syncParams();
            }
            null_value = std::max(m_logLike->value(), null_value);
            TsValues[srcNames[i]] = 2.*(logLike_value - null_value);
         } else {
// Null hypothesis has no remaining free parameters, so skip the fit
// (and hope that the model isn't empty).
            try {
               TsValues[srcNames[i]] = 2.*(logLike_value - m_logLike->value());
            } catch (std::exception & eObj) {
               m_formatter->err() << eObj.what() << std::endl;
            }
         }
         m_logLike->addSource(m_tsSrc);
         m_logLike->setFreeParamValues(fitParams);
      }
   }
   m_formatter->warn() << "!" << std::endl;
// Reset parameter values.
   m_logLike->setFreeParamValues(fitParams);
}

void likelihood::renormModel() {
   double freeNpred; 
   double totalNpred;
   npredValues(freeNpred, totalNpred);
   if (freeNpred <= 0) {
      return;
   }
   double deficit(observedCounts() - totalNpred);
   double renormFactor(1. + deficit/freeNpred);
   if (renormFactor < 1.) {
      renormFactor = 1.;
   }
   std::vector<std::string> srcNames;
   m_logLike->getSrcNames(srcNames);
   for (std::vector<std::string>::const_iterator srcName = srcNames.begin();
        srcName != srcNames.end(); ++srcName) {
      Source * src(m_logLike->getSource(*srcName));
      optimizers::Parameter & par(normPar(src));
      if (par.isFree() && isDiffuseOrNearby(src)) {
         double newValue(par.getValue()*renormFactor);
         par.setValue(newValue);
      }
   }   
}

double likelihood::observedCounts() {
   if (m_statistic == "BINNED") {
      double totalCounts(0);
      if (::getenv("USE_BINNED_LIKELIHOOD2")) {
         BinnedLikelihood2 * logLike = dynamic_cast<BinnedLikelihood2 *>(m_logLike);
         const std::vector<double> & counts(logLike->countsSpectrum());
         for (size_t i = 0; i < counts.size(); i++) {
            totalCounts += counts.at(i);
         }
      } else {
         BinnedLikelihood * logLike = dynamic_cast<BinnedLikelihood *>(m_logLike);
         const std::vector<double> & counts(logLike->countsSpectrum());
         for (size_t i = 0; i < counts.size(); i++) {
            totalCounts += counts.at(i);
         }
      }
      return totalCounts;
   }
   return m_helper->observation().eventCont().events().size();
}

void likelihood::npredValues(double & freeNpred, double & totalNpred) const {
   std::vector<std::string> srcNames;
   m_logLike->getSrcNames(srcNames);
   freeNpred = 0;
   totalNpred = 0;
   for (std::vector<std::string>::const_iterator srcName = srcNames.begin();
        srcName != srcNames.end(); ++srcName) {
      Source * src = m_logLike->getSource(*srcName);
      double npred(m_logLike->NpredValue(*srcName));
      totalNpred += npred;
      if (normPar(src).isFree() && isDiffuseOrNearby(src)) {
         freeNpred += npred;
      }
   }
}

optimizers::Parameter & likelihood::normPar(Source * src) const {
   const optimizers::Function & spectrum(src->spectrum());
   return const_cast<optimizers::Function &>(spectrum).normPar();
}

bool likelihood::isDiffuseOrNearby(Source * src) const {
   if (src->getType() == "Diffuse") {
      return true;
   } else if (::ptsrcSeparation(m_tsSrc, src) < m_maxdist) {
      return true;
   }
   return false;
}

bool likelihood::prompt(const std::string &query) {
   m_formatter->info(0) << query << std::endl;
   char answer[2];
   std::cin.getline(answer, 2);
   if (std::string(answer) == "y" || std::string(answer) == "") {
      return true;
   }
   return false;
}
