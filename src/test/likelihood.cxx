/**
 * @file likelihood.cxx
 * @brief Prototype standalone application for the Likelihood tool.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/test/likelihood.cxx,v 1.16 2003/12/05 05:39:23 jchiang Exp $
 */

#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <iostream>
#include <cstring>
#include <cmath>
#include <cassert>

#include "optimizers/FunctionFactory.h"
#include "optimizers/Lbfgs.h"
#include "optimizers/Minuit.h"
#include "optimizers/Drmngb.h"
#include "optimizers/Exception.h"

#include "latResponse/IrfsFactory.h"

#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/SourceModel.h" 
#include "Likelihood/Source.h"
#include "Likelihood/ScData.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/SpatialMap.h"
#include "Likelihood/LogLike.h"
#include "Likelihood/OptEM.h"
#include "Likelihood/RunParams.h"
#include "Likelihood/Exception.h"

using namespace Likelihood;

void print_fit_results(SourceModel &stat, const std::vector<double> &errors);
bool prompt(const std::string &query);

int main(int iargc, char* argv[]) {

#ifdef TRAP_FPE
   feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif

// Read in the command-line parameters using HOOPS
   std::string filename("likelihood.par");
   delete argv[0];
   argv[0] = strdup(filename.c_str());

   RunParams params(iargc, argv);

// Set the region-of-interest.
   std::string roiCutsFile;
   params.getParam("ROI_cuts_file", roiCutsFile);
   RoiCuts::setCuts(roiCutsFile);

// Read in the pointing information.
   std::string scFile;
   params.getParam("Spacecraft_file", scFile);
   long scHdu;
   params.getParam("Spacecraft_file_hdu", scHdu);
   std::vector<std::string> scFiles;
   RunParams::resolve_fits_files(scFile, scFiles);
   std::vector<std::string>::const_iterator scIt = scFiles.begin();
   for ( ; scIt != scFiles.end(); scIt++) {
      ScData::readData(*scIt, scHdu);
   }

// Read in the exposure map file.
   std::string exposureFile;
   params.getParam("Exposure_map_file", exposureFile);
   if (exposureFile != "none") {
      ExposureMap::readExposureFile(exposureFile);

   }

// Create the response functions.
   std::string responseFuncs;
   params.getParam("Response_functions", responseFuncs);
   latResponse::IrfsFactory irfsFactory;
   if (responseFuncs == "COMBINED_G25") {
      ResponseFunctions::addRespPtr(4, 
                                    irfsFactory.create("Glast25::Combined"));
   } else if (responseFuncs == "FRONT/BACK_G25") {
      ResponseFunctions::addRespPtr(2, irfsFactory.create("Glast25::Front"));
      ResponseFunctions::addRespPtr(3, irfsFactory.create("Glast25::Back"));
   } else if (responseFuncs == "TESTDC1") {
      ResponseFunctions::addRespPtr(1, irfsFactory.create("DC1::test"));
   } else if (responseFuncs == "FRONT") {
      ResponseFunctions::addRespPtr(5, irfsFactory.create("DC1::Front"));
   } else if (responseFuncs == "BACK") {
      ResponseFunctions::addRespPtr(6, irfsFactory.create("DC1::Back"));
   } else if (responseFuncs == "FRONT/BACK") {
      ResponseFunctions::addRespPtr(5, irfsFactory.create("DC1::Front"));
      ResponseFunctions::addRespPtr(6, irfsFactory.create("DC1::Back"));
   }

// Fill a FunctionFactory with Function object prototypes for source
// modeling.
   optimizers::FunctionFactory funcFactory;

// Add the prototypes for modeling spatial distributions.
   bool makeClone(false);
   funcFactory.addFunc("SkyDirFunction", new SkyDirFunction(), makeClone);
   funcFactory.addFunc("SpatialMap", new SpatialMap(), makeClone);

   LogLike * logLike = 0;
   bool useOptEM;
   params.getParam("Use_OptEM", useOptEM);
   if (useOptEM) {
      logLike = new OptEM();
   } else {
      logLike = new LogLike();
   }
   
// Read in the Event data.
   std::string eventFile;
   params.getParam("event_file", eventFile);
   long eventFileHdu;
   params.getParam("event_file_hdu", eventFileHdu);
   std::vector<std::string> eventFiles;
   RunParams::resolve_fits_files(eventFile, eventFiles);
   std::vector<std::string>::const_iterator evIt = eventFiles.begin();
   for ( ; evIt != eventFiles.end(); evIt++) {
      logLike->getEvents(*evIt, eventFileHdu);
   }

// Set the verbosity level and convergence tolerance.
   long verbose;
   params.getParam("fit_verbosity", verbose);
   double tol;
   params.getParam("fit_tolerance", tol);
   std::vector<double> errors;

// The fit loop.  If indicated, query the user at the end of each
// iteration whether the fit is to be performed again.  This allows
// the user to adjust the source model xml file by hand between
// iterations.

   bool queryLoop;
   params.getParam("query_for_refit", queryLoop);
   optimizers::Optimizer * myOpt = 0;
   do {
// Read in the Source model.
      std::string sourceModel;
      params.getParam("Source_model_file", sourceModel);
      if (logLike->getNumSrcs() == 0) {
// Read in the Source model for the first time.
         try {
            logLike->readXml(sourceModel, funcFactory);
         } catch (Likelihood::Exception &eObj) {
            std::cout << eObj.what();
            std::cout << "Check your source model file." << std::endl;
            assert(false);
         }
         logLike->computeEventResponses();
      } else {
// Re-read the Source model from the xml file, allowing only for 
// Parameter adjustments.
         logLike->reReadXml(sourceModel);
      }

// Do the fit.
      if (useOptEM) {
         dynamic_cast<OptEM *>(logLike)->findMin(verbose);
      } else {
// Select an optimizer.
         std::string optimizer;
         params.getParam("optimizer", optimizer);
         if (optimizer == "LBFGS") {
            myOpt = new optimizers::Lbfgs(*logLike);
         } else if (optimizer == "MINUIT") {
            myOpt = new optimizers::Minuit(*logLike);
         } else if (optimizer == "DRMNGB") {
            myOpt = new optimizers::Drmngb(*logLike);
         }
         try {
            myOpt->find_min(verbose, tol);
         } catch (optimizers::Exception &eObj) {
            std::cerr << eObj.what() << std::endl;
         }
// Evaluate the uncertainties.
         try {
            errors = myOpt->getUncertainty();
         } catch (optimizers::Exception &eObj) {
            std::cerr << "Error in computing uncertainties: \n"
                      << eObj.what() << "\n"
                      << "Bailing..."
                      << std::endl;
         } catch (...) {
            std::cerr << "Unexpected exception in computing uncertainties: \n"
                      << "Bailing..."
                      << std::endl;
         }
         if (optimizer == "DRMNGB") {
            int retCode =
               dynamic_cast<optimizers::Drmngb *>(myOpt)->getRetCode();
            std::cerr << "Drmngb return code: " << retCode;
         }
         delete myOpt;
      }
      print_fit_results(*logLike, errors);
      std::cout << std::endl;

// Write the model to the output xml file.
      std::string xmlFile;
      params.getParam("Source_model_output_file", xmlFile);
      std::string funcFileName("");
//      params.getParam("Function_models_file_name", funcFileName);

      if (xmlFile != "none") {
         std::cout << "Writing fitted model to " << xmlFile << std::endl;
         logLike->writeXml(xmlFile, funcFileName);
      }
   } while (queryLoop && prompt("Refit? [y] "));

// Write the model to a flux-style output file.
   std::string xml_fluxFile;
   params.getParam("flux_style_model_file", xml_fluxFile);
   if (xml_fluxFile != "none") {
      std::cout << "Writing flux-style xml model file to "
                << xml_fluxFile << std::endl;
      logLike->write_fluxXml(xml_fluxFile);
   }
   delete logLike;
}

void print_fit_results(SourceModel &stat, const std::vector<double> &errors) {
   std::vector<std::string> srcNames;
   stat.getSrcNames(srcNames);
   std::vector<optimizers::Parameter> parameters;

   std::vector<double>::const_iterator errIt = errors.begin();

   for (unsigned int i = 0; i < srcNames.size(); i++) {
      Source *src = stat.getSource(srcNames[i]);
      Source::FuncMap srcFuncs = src->getSrcFuncs();
      srcFuncs["Spectrum"]->getParams(parameters);
      std::cout << "\n" << srcNames[i] << ":\n";
      for (unsigned int i = 0; i < parameters.size(); i++) {
         std::cout << parameters[i].getName() << ": "
                   << parameters[i].getValue();
         if (parameters[i].isFree() && errIt != errors.end()) {
            std::cout << " +/- " << *errIt;
            errIt++;
         }
         std::cout << std::endl;
      }
      std::cout << "Npred: "
                << src->Npred() << std::endl;
   }
}

bool prompt(const std::string &query) {
   std::cout << query << std::endl;
   char answer[2];
   std::cin.getline(answer, 2);
   if (std::string(answer) == "y" || std::string(answer) == "") {
      return true;
   }
   return false;
}
