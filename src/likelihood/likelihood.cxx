/**
 * @file likelihood.cxx
 * @brief Prototype standalone application for the Likelihood tool.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/likelihood/likelihood.cxx,v 1.103 2006/02/23 01:54:53 jchiang Exp $
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

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "st_facilities/Util.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "optimizers/Drmngb.h"
#include "optimizers/Exception.h"
#include "optimizers/Lbfgs.h"
#include "optimizers/Minuit.h"
//#include "optimizers/newMinuit.h"
#ifdef HAVE_OPT_PP
#include "optimizers/OptPP.h"
#endif // HAVE_OPT_PP

#include "Likelihood/AppHelpers.h"
#include "Likelihood/BinnedLikelihood.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/ExposureCube.h"
#include "Likelihood/LogLike.h"
#include "Likelihood/MapShape.h"
#include "Likelihood/OptEM.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/Source.h"
#include "Likelihood/SourceMap.h"

#include "Verbosity.h"

namespace {
   class NormNames {
   public:
      static NormNames * instance() {
         if (s_instance == 0) {
            s_instance = new NormNames();
         }
         return s_instance;
      }

      const std::string & operator[](const std::string & name) const {
         std::map<std::string, std::string>::const_iterator it;
         it = m_normNames.find(name);
         return it->second;
      }
   private:
      static NormNames * s_instance;
      NormNames() {
         m_normNames["ConstantValue"] = "Value";
         m_normNames["BrokenPowerLaw"] = "Prefactor";
         m_normNames["BrokenPowerLaw2"] = "Integral";
         m_normNames["PowerLaw"] = "Prefactor";
         m_normNames["PowerLaw2"] = "Integral";
         m_normNames["Gaussian"] = "Prefactor";
         m_normNames["FileFunction"] = "Normalization";
         m_normNames["LogParabola"] = "norm";
      }
      std::map<std::string, std::string> m_normNames;
   };
   NormNames * NormNames::s_instance(0);

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
}

using namespace Likelihood;

/**
 * @class likelihood
 * @brief A class encapsulating the methods for performing a
 * Likelihood analysis in ballistic fashion.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/likelihood/likelihood.cxx,v 1.103 2006/02/23 01:54:53 jchiang Exp $
 */

class likelihood : public st_app::StApp {
public:
   likelihood();
   virtual ~likelihood() throw() {
      try {
         delete m_logLike;
         delete m_opt;
         delete m_dataMap;
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
   std::string m_statistic;

   std::clock_t m_cpuStart;

   void promptForParameters();
   void createStatistic();
   void readEventData();
   void readSourceModel();
   void selectOptimizer(std::string optimizer="");
   void writeSourceXml();
   void writeFluxXml();
   void writeCountsSpectra();
   void writeCountsMap();
   void printFitResults(const std::vector<double> &errors);
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

   static std::string s_cvs_id;
};

st_app::StAppFactory<likelihood> myAppFactory("gtlikelihood");

std::string likelihood::s_cvs_id("$Name:  $");

void likelihood::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

likelihood::likelihood() 
   : st_app::StApp(), m_helper(0), 
     m_pars(st_app::StApp::getParGroup("gtlikelihood")),
     m_logLike(0), m_opt(0), m_dataMap(0), m_cpuStart(std::clock()),
     m_tsSrc(0), m_maxdist(20.) {
   setVersion(s_cvs_id);
   m_pars.setSwitch("statistic");
   m_pars.setCase("statistic", "BINNED", "counts_map_file");
   m_pars.setCase("statistic", "BINNED", "binned_exposure_map");
   m_pars.setCase("statistic", "BINNED", "apply_psf_corrections");
   m_pars.setCase("statistic", "UNBINNED", "evfile");
   m_pars.setCase("statistic", "UNBINNED", "evtable");
   m_pars.setCase("statistic", "UNBINNED", "scfile");
   m_pars.setCase("statistic", "UNBINNED", "sctable");
   m_pars.setCase("statistic", "UNBINNED", "exposure_map_file");
}

void likelihood::run() {
   promptForParameters();
   Verbosity::instance(m_pars["chatter"]);
   m_helper = new AppHelpers(&m_pars);
   std::string expcube_file = m_pars["exposure_cube_file"];
   if (expcube_file != "none" && expcube_file != "") {
      ExposureCube & expCube = 
         const_cast<ExposureCube &>(m_helper->observation().expCube());
      expCube.readExposureCube(expcube_file);
   }
   bool useEdisp = m_pars["use_energy_dispersion"];
   ResponseFunctions & respFuncs = 
      const_cast<ResponseFunctions &>(m_helper->observation().respFuncs());
   respFuncs.setEdispFlag(useEdisp);
   if (m_statistic == "BINNED") {
      m_helper->setRoi(m_pars["counts_map_file"], "", false);
   } else {
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
   }
   createStatistic();

// Set the verbosity level and convergence tolerance.
   long verbose = m_pars["chatter"];
// ST proclaimed nominal verbosity level is 2, but optimizers expect 1, so
// we subtract 1.
   if (verbose > 1) verbose--;
   double tol = m_pars["fit_tolerance"];
   std::vector<double> errors;

// The fit loop.  If indicated, query the user at the end of each
// iteration whether the fit is to be performed again.  This allows
// the user to adjust the source model xml file by hand between
// iterations.
   bool queryLoop = m_pars["query_for_refit"];
   do {
      readSourceModel();
// Do the fit.
/// @todo Allow the optimizer to be re-selected here by the user.    
      selectOptimizer();
      try {
         m_opt->find_min(verbose, tol);
         try {
            errors = m_opt->getUncertainty();
            setErrors(errors);
         } catch (optimizers::Exception & eObj) {
            std::cerr << "Exception encountered while estimating errors:\n";
            std::cerr << eObj.what() << std::endl;
//             throw;
         }
      } catch (optimizers::Exception & eObj) {
         std::cerr << "Exception encountered while minimizing "
                   << "objective function:\n";
         std::cerr << eObj.what() << std::endl;
//          throw;
      }
      printFitResults(errors);
      writeSourceXml();
   } while (queryLoop && prompt("Refit? [y] "));
   writeFluxXml();
   if (m_statistic != "BINNED") {
      writeCountsSpectra();
   }
//   writeCountsMap();
   if (print_output()) {
      std::cout << "Elapsed CPU time: " << cputime() << std::endl;
   }
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
      m_pars.Prompt("counts_map_file");
      m_pars.Prompt("binned_exposure_map");
   } else {
      m_pars.Prompt("scfile");
      m_pars.Prompt("evfile");
      m_pars.Prompt("exposure_map_file");
   }
   m_pars.Prompt("exposure_cube_file");
   m_pars.Prompt("source_model_file");
   m_pars.Prompt("source_model_output_file");
   AppHelpers::checkOutputFile(m_pars["clobber"], 
                               m_pars["source_model_output_file"]);
   m_pars.Prompt("flux_style_model_file");
   AppHelpers::checkOutputFile(m_pars["clobber"], 
                               m_pars["flux_style_model_file"]);
   m_pars.Prompt("rspfunc");
   m_pars.Prompt("use_energy_dispersion");
   m_pars.Prompt("optimizer");
   m_pars.Prompt("query_for_refit");

   m_pars.Save();
}

void likelihood::createStatistic() {
   if (m_statistic == "BINNED") {
      if (!m_helper->observation().expCube().haveFile()) {
         throw std::runtime_error
            ("An exposure cube file is required for binned analysis. "
             "Please specify an exposure cube file.");
      }
      std::string countsMapFile = m_pars["counts_map_file"];
      st_facilities::Util::file_ok(countsMapFile);
      m_dataMap = new CountsMap(countsMapFile);
      bool apply_psf_corrections(false);
      try {
         apply_psf_corrections = m_pars["apply_psf_corrections"];
      } catch (...) {
         // assume parameter does not exist, so use default value.
      }
      m_logLike = new BinnedLikelihood(*m_dataMap, m_helper->observation(),
                                       countsMapFile, apply_psf_corrections);
      std::string binnedMap = m_pars["binned_exposure_map"];
      if (binnedMap != "none" && binnedMap != "") {
         SourceMap::setBinnedExposure(binnedMap);
      }
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
      m_helper->observation().eventCont().getEvents(*evIt);
   }
}

void likelihood::readSourceModel() {
   std::string sourceModel = m_pars["source_model_file"];
   bool requireExposure(true);
   if (m_statistic == "BINNED") {
      requireExposure = false;
   }
   if (m_logLike->getNumSrcs() == 0) {
// Read in the Source model for the first time.
      st_facilities::Util::file_ok(sourceModel);
      m_logLike->readXml(sourceModel, m_helper->funcFactory(), 
                         requireExposure);
      if (m_statistic != "BINNED") {
         m_logLike->computeEventResponses();
      } else {
//         dynamic_cast<BinnedLikelihood *>(m_logLike)->saveSourceMaps();
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
   if (optimizer == "LBFGS") {
      m_opt = new optimizers::Lbfgs(*m_logLike);
   } else if (optimizer == "MINUIT") {
      m_opt = new optimizers::Minuit(*m_logLike);
//    } else if (optimizer == "NEWMINUIT") {
//       m_opt = new optimizers::newMinuit(*m_logLike);
   } else if (optimizer == "DRMNGB") {
      m_opt = new optimizers::Drmngb(*m_logLike);
#ifdef HAVE_OPT_PP
   } else if (optimizer == "OPTPP") {
      m_opt = new optimizers::OptPP(*m_logLike);
#endif // HAVE_OPT_PP
   }
   if (m_opt == 0) {
      throw std::invalid_argument("Invalid optimizer choice: " + optimizer);
   }
}

void likelihood::writeSourceXml() {
   std::string xmlFile = m_pars["source_model_output_file"];
   std::string funcFileName("");
   if (xmlFile != "none") {
      if (print_output()) {
         std::cout << "Writing fitted model to " << xmlFile << std::endl;
      }
      m_logLike->writeXml(xmlFile, funcFileName);
   }
}

void likelihood::writeFluxXml() {
   std::string xml_fluxFile = m_pars["flux_style_model_file"];
   if (xml_fluxFile != "none") {
      if (print_output()) {
         std::cout << "Writing flux-style xml model file to "
                   << xml_fluxFile << std::endl;
      }
      m_logLike->write_fluxXml(xml_fluxFile);
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
   const EventData myData(m_helper->observation().eventCont().events());
   std::vector<double> energies;
   const RoiCuts & roiCuts = m_helper->observation().roiCuts();
   double emin = roiCuts.getEnergyCuts().first;
   double emax = roiCuts.getEnergyCuts().second;
   int nee(20);
   double estep = log(emax/emin)/(nee-1);
   for (int k = 0; k < nee; k++) {
      energies.push_back(emin*exp(estep*k));
   }
   std::vector<double> evals;
   std::vector<std::string> srcNames;
   m_logLike->getSrcNames(srcNames);
   std::vector< std::vector<double> > npred(srcNames.size());
   std::ofstream outputFile("counts.dat");
   for (int k = 0; k < nee - 1; k++) {
      bool writeLine(true);
      std::ostringstream line;
      line << sqrt(energies[k]*energies[k+1]) << "   "
           << myData.nobs(energies[k], energies[k+1]) << "  ";
      for (unsigned int i = 0; i < srcNames.size(); i++) {
         Source * src = m_logLike->getSource(srcNames[i]);
         double Npred;
         try {
            Npred = src->Npred(energies[k], energies[k+1]);
            if (i==0) evals.push_back(log10(sqrt(energies[k]*energies[k+1])));
            line << Npred << "  ";
         } catch (std::out_of_range &) {
            writeLine = false;
         }
      }
      if (writeLine) {
         outputFile << line.str() << std::endl;
      }
   }
   outputFile.close();
}

void likelihood::writeCountsMap() {
// If there is no valid exposure_cube_file, do nothing and return.
   std::string expcube_file = m_pars["exposure_cube_file"];
   if (expcube_file == "none") {
      return;
   }
   ExposureCube & expCube = 
      const_cast<ExposureCube &>(m_helper->observation().expCube());
   expCube.readExposureCube(expcube_file);
   m_dataMap->writeOutput("likelihood", "data_map.fits");
   try {
      CountsMap * modelMap;
      if (m_statistic == "BINNED") {
         modelMap = m_logLike->createCountsMap();
      } else {
         modelMap = m_logLike->createCountsMap(*m_dataMap);
      }
      modelMap->writeOutput("likelihood", "model_map.fits");
      delete modelMap;
   } catch (std::exception & eObj) {
      std::cout << eObj.what() << std::endl;
   }
}

void likelihood::printFitResults(const std::vector<double> &errors) {
   std::vector<std::string> srcNames;
   m_logLike->getSrcNames(srcNames);
   
   std::map<std::string, double> TsValues;
   std::map<std::string, double> RoiDist;

   computeTsValues(srcNames, TsValues, RoiDist);

   std::vector<optimizers::Parameter> parameters;
   std::vector<double>::const_iterator errIt = errors.begin();

   std::ofstream resultsFile("results.dat");
   resultsFile << "{";

   double totalNpred(0);
   for (unsigned int i = 0; i < srcNames.size(); i++) {
      Source * src = m_logLike->getSource(srcNames[i]);
      Source::FuncMap srcFuncs = src->getSrcFuncs();
      srcFuncs["Spectrum"]->getParams(parameters);
      if (print_output()) {
         std::cout << "\n" << srcNames[i] << ":\n";
      }
      resultsFile << "'" << srcNames[i] << "': {";
      for (unsigned int j = 0; j < parameters.size(); j++) {
         if (print_output()) {
            std::cout << parameters[j].getName() << ": "
                      << parameters[j].getValue();
         }
         resultsFile << "'" << parameters[j].getName() << "': "
                     << "'" << parameters[j].getValue();
         if (parameters[j].isFree() && errIt != errors.end()) {
            if (print_output()) {
               std::cout << " +/- " << *errIt;
            }
            resultsFile << " +/- " << *errIt << "',\n";
            errIt++;
         } else {
            resultsFile << "',\n";
         }
         if (print_output()) {
            std::cout << std::endl;
         }
      }
      if (m_statistic != "BINNED") {
         double npred(src->Npred());
         if (print_output()) {
            std::cout << "Npred: "
                      << npred << std::endl;
         }
         resultsFile << "'Npred': '" << npred << "',\n";
         totalNpred += npred;
      }
      if (RoiDist.count(srcNames[i])) {
         if (print_output()) {      
            std::cout << "ROI distance: "
                      << RoiDist[srcNames[i]] << std::endl;
         }
         resultsFile << "'ROI distance': '" << RoiDist[srcNames[i]] << "',\n";
      }
      if (TsValues.count(srcNames[i])) {
         if (print_output()) {
            std::cout << "TS value: "
                      << TsValues[srcNames[i]] << std::endl;
         }
         resultsFile << "'TS value': '" << TsValues[srcNames[i]] << "',\n";
      }
      resultsFile << "},\n";
   }
   resultsFile << "}" << std::endl;
   if (print_output()) {
      std::cout << "\nTotal number of observed counts: "
                << observedCounts() << std::endl;
      std::cout << "Total number of model events: ";
      if (m_statistic == "BINNED") {
         std::cout << dynamic_cast<BinnedLikelihood *>(m_logLike)->npred();
      } else {
         std::cout << totalNpred;
      }
      std::cout << std::endl;
   }
   resultsFile.close();
   if (print_output()) {
      std::cout << "\n-log(Likelihood): "
                << std::setprecision(10)
                << -m_logLike->value()
                << "\n" << std::endl;
   }
   delete m_opt;
   m_opt = 0;
}

void likelihood::computeTsValues(const std::vector<std::string> & srcNames,
                                 std::map<std::string, double> & TsValues,
                                 std::map<std::string, double> & RoiDist) {
// Save current set of parameters.
   std::vector<double> fitParams;
   m_logLike->getFreeParamValues(fitParams);
   double logLike_value = m_logLike->value();

   int verbose(0);
   double tol(1e-4);
   if (print_output()) {
      std::cerr << "Computing TS values for each source ("
                << srcNames.size() << " total)\n";
   }
   astro::SkyDir roiCenter 
      = m_helper->observation().roiCuts().extractionRegion().center();
   for (unsigned int i = 0; i < srcNames.size(); i++) {
      if (print_output()) {
         std::cerr << ".";
      }
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
            if (m_pars["find_Ts_mins"]) {
               try {
                  m_opt->find_min(verbose, tol);
               } catch (std::exception & eObj) {
                  std::cout << eObj.what() << std::endl;
               }
            } else {
               if (m_statistic != "BINNED") {
//                  renormModel();
               }
            }
            null_value = std::max(m_logLike->value(), null_value);
            TsValues[srcNames[i]] = 2.*(logLike_value - null_value);
         } else {
// Null hypothesis has no remaining free parameters, so skip the fit
// (and hope that the model isn't empty).
            try {
               TsValues[srcNames[i]] = 2.*(logLike_value - m_logLike->value());
            } catch (std::exception & eObj) {
               std::cout << eObj.what() << std::endl;
            }
         }
         m_logLike->addSource(m_tsSrc);
         m_logLike->setFreeParamValues(fitParams);
      }
   }
   if (print_output()) {
      std::cerr << "!" << std::endl;
   }
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
      BinnedLikelihood * logLike = dynamic_cast<BinnedLikelihood *>(m_logLike);
      const std::vector<double> & counts(logLike->countsSpectrum());
      double totalCounts(0);
      for (size_t i = 0; i < counts.size(); i++) {
         totalCounts += counts.at(i);
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
      double npred(src->Npred());
      totalNpred += npred;
      if (normPar(src).isFree() && isDiffuseOrNearby(src)) {
         freeNpred += npred;
      }
   }
}

optimizers::Parameter & likelihood::normPar(Source * src) const {
   ::NormNames & normNames(*::NormNames::instance());
   const optimizers::Function & spectrum(src->spectrum());
   const std::string & parname(normNames[spectrum.genericName()]);
   return const_cast<optimizers::Function &>(spectrum).parameter(parname);
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
   std::cout << query << std::endl;
   char answer[2];
   std::cin.getline(answer, 2);
   if (std::string(answer) == "y" || std::string(answer) == "") {
      return true;
   }
   return false;
}
