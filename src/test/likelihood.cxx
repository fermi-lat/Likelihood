/**
 * @file likelihood.cxx
 * @brief Prototype standalone application for the Likelihood tool.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/test/likelihood.cxx,v 1.25 2004/03/17 19:09:54 jchiang Exp $
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

#include "hoops/hoops.h"
#include "hoops/hoops_prompt_group.h"

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
#include "Likelihood/Exception.h"

using namespace Likelihood;
using latResponse::irfsFactory;

void print_fit_results(LogLike &logLike, const std::vector<double> &errors);
bool prompt(const std::string &query);

namespace {
   bool fileExists(const std::string &filename) {
      std::ifstream file(filename.c_str());
      return file.is_open();
   }

   void file_ok(std::string filename) {
      facilities::Util::expandEnvVar(&filename);
      if (fileExists(filename)) {
         return;
      } else {
         std::cout << "likelihood::main:\n"
                   << "File not found: " << filename
                   << std::endl;
//         assert(::fileExists(filename));
         exit(-1);
      }
   }

   void readLines(std::string inputFile, 
                  std::vector<std::string> &lines) {
      
      facilities::Util::expandEnvVar(&inputFile);
      
      std::ifstream file(inputFile.c_str());
      lines.clear();
      std::string line;
      while (std::getline(file, line, '\n')) {
         if (line != "" && line != " ") { //skip (most) blank lines
            lines.push_back(line);
         }
      }
   }

   void resolve_fits_files(std::string filename, 
                           std::vector<std::string> &files) {
      
      facilities::Util::expandEnvVar(&filename);
      files.clear();
      
// Read the first line of the file and see if the first 6 characters
// are "SIMPLE".  If so, then we assume it's a FITS file.
      std::ifstream file(filename.c_str());
      std::string firstLine;
      std::getline(file, firstLine, '\n');
      if (firstLine.find("SIMPLE") == 0) {
// This is a FITS file. Return that as the sole element in the files
// vector.
         files.push_back(filename);
         return;
      } else {
// filename contains a list of fits files.
         readLines(filename, files);
         return;
      }
   }

}

int main(int iargc, char* argv[]) {

#ifdef TRAP_FPE
   feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif

// Read in the command-line parameters using HOOPS
   strcpy(argv[0], "likelihood");

   try {
      hoops::ParPromptGroup pars(iargc, argv);
      pars.Prompt();

// Set the region-of-interest.
      std::string roiCutsFile = pars["ROI_cuts_file"];
      ::file_ok(roiCutsFile);
      RoiCuts::setCuts(roiCutsFile);

// Read in the pointing information.
      std::string scFile = pars["Spacecraft_file"];
      ::file_ok(scFile);
      long scHdu = pars["Spacecraft_file_hdu"];
      std::vector<std::string> scFiles;
      ::resolve_fits_files(scFile, scFiles);
      std::vector<std::string>::const_iterator scIt = scFiles.begin();
      for ( ; scIt != scFiles.end(); scIt++) {
         ::file_ok(*scIt);
         ScData::readData(*scIt, scHdu);
      }

// Read in the exposure map file.
      std::string exposureFile = pars["Exposure_map_file"];
      if (exposureFile != "none") {
         ::file_ok(exposureFile);
         ExposureMap::readExposureFile(exposureFile);
      }

// Create the response functions.
      std::string responseFuncs = pars["Response_functions"];
      std::map< std::string, std::vector<std::string> > responseIds;
      responseIds["FRONT"].push_back("DC1::Front");
      responseIds["BACK"].push_back("DC1::Back");
      responseIds["FRONT/BACK"].push_back("DC1::Front");
      responseIds["FRONT/BACK"].push_back("DC1::Back");
      responseIds["GLAST25"].push_back("Glast25::Front");
      responseIds["GLAST25"].push_back("Glast25::Back");

      if (responseIds.count(responseFuncs)) {
         std::vector<std::string> &resps = responseIds[responseFuncs];
         for (unsigned int i = 0; i < resps.size(); i++) {
            ResponseFunctions::addRespPtr(i, irfsFactory().create(resps[i]));
         }
      } else {
         std::cerr << "Invalid response function choice: "
                   << responseFuncs << std::endl;
         exit(-1);
      }

// Fill a FunctionFactory with Function object prototypes for source
// modeling.
      optimizers::FunctionFactory funcFactory;

// Add the prototypes for modeling spatial distributions.
      bool makeClone(false);
      funcFactory.addFunc("SkyDirFunction", new SkyDirFunction(), makeClone);
      funcFactory.addFunc("SpatialMap", new SpatialMap(), makeClone);
      
      LogLike * logLike = 0;
      bool useOptEM = pars["Use_OptEM"];
      if (useOptEM) {
         logLike = new OptEM();
      } else {
         logLike = new LogLike();
      }
   
// Read in the Event data.
      std::string eventFile = pars["event_file"];
      long eventFileHdu = pars["event_file_hdu"];
      std::vector<std::string> eventFiles;
      ::file_ok(eventFile);
      ::resolve_fits_files(eventFile, eventFiles);
      std::vector<std::string>::const_iterator evIt = eventFiles.begin();
      for ( ; evIt != eventFiles.end(); evIt++) {
         ::file_ok(*evIt);
         logLike->getEvents(*evIt, eventFileHdu);
      }

// Set the verbosity level and convergence tolerance.
      long verbose = pars["fit_verbosity"];
      double tol = pars["fit_tolerance"];
      std::vector<double> errors;

// The fit loop.  If indicated, query the user at the end of each
// iteration whether the fit is to be performed again.  This allows
// the user to adjust the source model xml file by hand between
// iterations.

      bool queryLoop = pars["query_for_refit"];
      optimizers::Optimizer * myOpt = 0;
      do {
// Read in the Source model.
         std::string sourceModel = pars["Source_model_file"];
         if (logLike->getNumSrcs() == 0) {
// Read in the Source model for the first time.
            try {
               ::file_ok(sourceModel);
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
            ::file_ok(sourceModel);
            logLike->reReadXml(sourceModel);
         }

// Do the fit.
         if (useOptEM) {
            try {
               dynamic_cast<OptEM *>(logLike)->findMin(verbose);
            } catch (optimizers::Exception &eObj) {
               std::cerr << eObj.what() << std::endl;
            }
         } else {
// Select an optimizer.
            std::string optimizer = pars["optimizer"];
            try {
               if (optimizer == "LBFGS") {
                  myOpt = new optimizers::Lbfgs(*logLike);
               } else if (optimizer == "MINUIT") {
                  myOpt = new optimizers::Minuit(*logLike);
               } else if (optimizer == "DRMNGB") {
                  myOpt = new optimizers::Drmngb(*logLike);
               }
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
               std::cerr << "Unexpected exception in "
                         << "computing uncertainties: \n"
                         << "Bailing..."
                         << std::endl;
            }
            if (optimizer == "DRMNGB") {
               int retCode =
                  dynamic_cast<optimizers::Drmngb *>(myOpt)->getRetCode();
               std::cerr << "Drmngb return code: " << retCode;
            }
            delete myOpt;
         } // useOptEM
         print_fit_results(*logLike, errors);
         std::cout << std::endl 
                   << "-log(Likelihood): " << -logLike->value()
                   << std::endl;
         std::cout << std::endl;
         
// Write the model to the output xml file.
         std::string xmlFile = pars["Source_model_output_file"];
         std::string funcFileName("");
//          std::string funcFileName = pars["Function_models_file_name"];

         if (xmlFile != "none") {
            std::cout << "Writing fitted model to " << xmlFile << std::endl;
            logLike->writeXml(xmlFile, funcFileName);
         }
      } while (queryLoop && prompt("Refit? [y] "));

// Write the model to a flux-style output file.
      std::string xml_fluxFile = pars["flux_style_model_file"];
      if (xml_fluxFile != "none") {
         std::cout << "Writing flux-style xml model file to "
                   << xml_fluxFile << std::endl;
         logLike->write_fluxXml(xml_fluxFile);
      }
      delete logLike;
   } catch (std::exception &eObj) {
      std::cout << eObj.what() << std::endl;
   }
}

void print_fit_results(LogLike &logLike, const std::vector<double> &errors) {
   std::vector<std::string> srcNames;
   logLike.getSrcNames(srcNames);

// Compute TS for each source.
   std::map<std::string, double> TsValues;
   int verbose(0);
   double tol(1e-4);
   double logLike_value = logLike.value();
   for (unsigned int i = 0; i < srcNames.size(); i++) {
      if (srcNames[i].find("Diffuse") == std::string::npos) {
         Source * src = logLike.deleteSource(srcNames[i]);
         if (logLike.getNumFreeParams() > 0) {
// Don't fit if there are no free parameters remaining.
            optimizers::Drmngb opt(logLike);
            opt.find_min(verbose, tol);
            TsValues[srcNames[i]] = 2.*(logLike_value - logLike.value());
         } else {
// // Not sure this is correct in the case where the model for the null
// // hypothesis is truly empty.
//             TsValues[srcNames[i]] = 2.*logLike_value;
// A better default value?
            TsValues[srcNames[i]] = 0.;
         }            
         logLike.addSource(src);
      }
   }

// Restore parameters to their previously fitted values.
   optimizers::Drmngb opt(logLike);
   opt.find_min(verbose, tol);

   std::vector<optimizers::Parameter> parameters;
   std::vector<double>::const_iterator errIt = errors.begin();

   for (unsigned int i = 0; i < srcNames.size(); i++) {
      Source *src = logLike.getSource(srcNames[i]);
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
