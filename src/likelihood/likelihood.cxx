/**
 * @file likelihood.cxx
 * @brief Prototype standalone application for the Likelihood tool.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/likelihood/likelihood.cxx,v 1.2 2004/04/06 03:32:26 jchiang Exp $
 */

#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <cmath>
#include <cstring>
#include <cassert>
#include <fstream>
#include <iostream>

#include "facilities/Util.h"

#include "st_app/IApp.h"
#include "hoops/hoops_prompt_group.h"

#include "optimizers/FunctionFactory.h"
#include "optimizers/Lbfgs.h"
#include "optimizers/Minuit.h"
#include "optimizers/Drmngb.h"
#include "optimizers/Exception.h"

#include "latResponse/IrfsFactory.h"

#include "Likelihood/AppBase.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/SourceModel.h" 
#include "Likelihood/Source.h"
#include "Likelihood/ScData.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/SpatialMap.h"
#include "Likelihood/LogLike.h"
#include "Likelihood/OptEM.h"
#include "Likelihood/Exception.h"
#include "Likelihood/Util.h"

using namespace Likelihood;
using latResponse::irfsFactory;

class likelihood : public AppBase {

public:

   likelihood() : AppBase("likelihood"), m_logLike(0),
                  m_useOptEM(false), m_opt(0) {}

   virtual ~likelihood() throw() {}

   virtual void run();

private:

   LogLike * m_logLike;
   bool m_useOptEM;
   optimizers::Optimizer * m_opt;

   virtual void tearDown();
   void createStatistic();
   void readEventData();
   void readSourceModel();
   void selectOptimizer();
   void writeSourceXml();
   void writeFluxXml();
   void printFitResults(const std::vector<double> &errors);
   bool prompt(const std::string &query);

};

st_app::IApp * my_application = new likelihood();

void likelihood::run() {
   setUp();
   readExposureMap();
   createStatistic();
   readEventData();

   hoops::IParGroup & pars = hoopsGetParGroup();

// Set the verbosity level and convergence tolerance.
   long verbose = pars["fit_verbosity"];
   double tol = pars["fit_tolerance"];
   std::vector<double> errors;

// The fit loop.  If indicated, query the user at the end of each
// iteration whether the fit is to be performed again.  This allows
// the user to adjust the source model xml file by hand between
// iterations.
   bool queryLoop = pars["query_for_refit"];
   do {
      readSourceModel();
// Do the fit.
      if (m_useOptEM) {
         dynamic_cast<OptEM *>(m_logLike)->findMin(verbose);
      } else {
// @todo Allow the optimizer to be re-selected here.         
         selectOptimizer();
         m_opt->find_min(verbose, tol);
         errors = m_opt->getUncertainty();
         delete m_opt;
         m_opt = 0;
      }
      printFitResults(errors);
      writeSourceXml();
   } while (queryLoop && prompt("Refit? [y] "));
   writeFluxXml();
   tearDown();
}

void likelihood::tearDown() {
   delete m_logLike;
}

void likelihood::createStatistic() {
   hoops::IParGroup & pars = hoopsGetParGroup();
   m_useOptEM = pars["Use_OptEM"];
   if (m_useOptEM) {
      m_logLike = new OptEM();
   } else {
      m_logLike = new LogLike();
   }
}

void likelihood::readEventData() {
   hoops::IParGroup & pars = hoopsGetParGroup();
   std::string eventFile = pars["event_file"];
   long eventFileHdu = pars["event_file_hdu"];
   std::vector<std::string> eventFiles;
   Util::file_ok(eventFile);
   Util::resolve_fits_files(eventFile, eventFiles);
   std::vector<std::string>::const_iterator evIt = eventFiles.begin();
   for ( ; evIt != eventFiles.end(); evIt++) {
      Util::file_ok(*evIt);
      m_logLike->getEvents(*evIt, eventFileHdu);
   }
}

void likelihood::readSourceModel() {
   hoops::IParGroup & pars = hoopsGetParGroup();
   std::string sourceModel = pars["Source_model_file"];
   if (m_logLike->getNumSrcs() == 0) {
// Read in the Source model for the first time.
      Util::file_ok(sourceModel);
      m_logLike->readXml(sourceModel, m_funcFactory);
      m_logLike->computeEventResponses();
   } else {
// Re-read the Source model from the xml file, allowing only for 
// Parameter adjustments.
      Util::file_ok(sourceModel);
      m_logLike->reReadXml(sourceModel);
   }
}

void likelihood::selectOptimizer() {
   hoops::IParGroup & pars = hoopsGetParGroup();
   std::string optimizer = pars["optimizer"];
   if (optimizer == "LBFGS") {
      m_opt = new optimizers::Lbfgs(*m_logLike);
   } else if (optimizer == "MINUIT") {
      m_opt = new optimizers::Minuit(*m_logLike);
   } else if (optimizer == "DRMNGB") {
      m_opt = new optimizers::Drmngb(*m_logLike);
   }
   if (m_opt == 0) {
      throw std::invalid_argument("Invalid optimizer choice: " + optimizer);
   }
}

void likelihood::writeSourceXml() {
   hoops::IParGroup & pars = hoopsGetParGroup();
   std::string xmlFile = pars["Source_model_output_file"];
   std::string funcFileName("");
   if (xmlFile != "none") {
      std::cout << "Writing fitted model to " << xmlFile << std::endl;
      m_logLike->writeXml(xmlFile, funcFileName);
   }
}

void likelihood::writeFluxXml() {
   hoops::IParGroup & pars = hoopsGetParGroup();
   std::string xml_fluxFile = pars["flux_style_model_file"];
   if (xml_fluxFile != "none") {
      std::cout << "Writing flux-style xml model file to "
                << xml_fluxFile << std::endl;
      m_logLike->write_fluxXml(xml_fluxFile);
   }
}

void likelihood::printFitResults(const std::vector<double> &errors) {
   std::vector<std::string> srcNames;
   m_logLike->getSrcNames(srcNames);

// Compute TS for each source.
   std::map<std::string, double> TsValues;
   int verbose(0);
   double tol(1e-4);
   double logLike_value = m_logLike->value();
   for (unsigned int i = 0; i < srcNames.size(); i++) {
      if (srcNames[i].find("Diffuse") == std::string::npos) {
         Source * src = m_logLike->deleteSource(srcNames[i]);
         if (m_logLike->getNumFreeParams() > 0) {
// Don't fit if there are no free parameters remaining.
            optimizers::Drmngb opt(*m_logLike);
            opt.find_min(verbose, tol);
            TsValues[srcNames[i]] = 2.*(logLike_value - m_logLike->value());
         } else {
// // Not sure this is correct in the case where the model for the null
// // hypothesis is truly empty.
//             TsValues[srcNames[i]] = 2.*logLike_value;
// A better default value?
            TsValues[srcNames[i]] = 0.;
         }            
         m_logLike->addSource(src);
      }
   }

// Restore parameters to their previously fitted values.
   optimizers::Drmngb opt(*m_logLike);
   opt.find_min(verbose, tol);

   std::vector<optimizers::Parameter> parameters;
   std::vector<double>::const_iterator errIt = errors.begin();

   for (unsigned int i = 0; i < srcNames.size(); i++) {
      Source *src = m_logLike->getSource(srcNames[i]);
      Source::FuncMap srcFuncs = src->getSrcFuncs();
      srcFuncs["Spectrum"]->getParams(parameters);
      std::cout << "\n" << srcNames[i] << ":\n";
      for (unsigned int j = 0; j < parameters.size(); j++) {
         std::cout << parameters[j].getName() << ": "
                   << parameters[j].getValue();
         if (parameters[j].isFree() && errIt != errors.end()) {
            std::cout << " +/- " << *errIt;
            errIt++;
         }
         std::cout << std::endl;
      }
      std::cout << "Npred: "
                << src->Npred() << std::endl;
      if (TsValues.count(srcNames[i])) {
         std::cout << "TS value: "
                   << TsValues[srcNames[i]] << std::endl;
      }
   }
   std::cout << std::endl 
             << "-log(Likelihood): " << -m_logLike->value()
             << std::endl;
   std::cout << std::endl;
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
