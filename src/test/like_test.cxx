/**
 * @file like_test.cxx
 * @brief Prototype standalone application for the Likelihood tool.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/test/like_test.cxx,v 1.1 2003/11/08 02:55:56 jchiang Exp $
 */

#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>

#include "facilities/Util.h"

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
#include "Likelihood/ConstantValue.h"
#include "Likelihood/LogLike.h"
#include "Likelihood/OptEM.h"
#include "Likelihood/RunParams.h"
#include "BrokenPowerLaw.h"
#include "PowerLaw.h"
#include "Gaussian.h"
#include "AbsEdge.h"

using namespace Likelihood;

void print_fit_results(SourceModel &stat, const std::vector<double> &errors);
void resolve_fits_files(std::string filename, std::vector<std::string> &files);
void readLines(std::string inputFile, std::vector<std::string> &lines);

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
   std::string roiCutsFile = params.string_par("ROI_cuts_file");
   RoiCuts::setCuts(roiCutsFile);

// Read in the pointing information.
   std::string scFile = params.string_par("Spacecraft_file");
   int scHdu = static_cast<int>(params.long_par("Spacecraft_file_hdu"));
   std::vector<std::string> scFiles;
   resolve_fits_files(scFile, scFiles);
   std::vector<std::string>::const_iterator scIt = scFiles.begin();
   for ( ; scIt != scFiles.end(); scIt++) {
      ScData::readData(*scIt, scHdu);
   }

// Read in the exposure map file.
   std::string exposureFile = params.string_par("Exposure_map_file");
   if (exposureFile != "none") {
      ExposureMap::readExposureFile(exposureFile);
   }

// Create the response functions.
   std::string responseFuncs = params.string_par("Response_functions");
   latResponse::IrfsFactory irfsFactory;
   if (responseFuncs == "COMBINED") {
      ResponseFunctions::addRespPtr(4, 
                                    irfsFactory.create("Glast25::Combined"));
   } else if (responseFuncs == "FRONT/BACK") {
      ResponseFunctions::addRespPtr(2, irfsFactory.create("Glast25::Front"));
      ResponseFunctions::addRespPtr(3, irfsFactory.create("Glast25::Back"));
   }

// Fill a FunctionFactory with Function object prototypes for source
// modeling.
   optimizers::FunctionFactory funcFactory;

// Add the standard prototypes for modeling spectra,
   bool makeClone(false);
   funcFactory.addFunc("PowerLaw", new PowerLaw(), makeClone);
   funcFactory.addFunc("BrokenPowerLaw", new BrokenPowerLaw(), makeClone);
   funcFactory.addFunc("Gaussian", new Gaussian(), makeClone);
   funcFactory.addFunc("AbsEdge", new AbsEdge(), makeClone);

// and the prototypes for modeling spatial distributions.
   funcFactory.addFunc("SkyDirFunction", new SkyDirFunction(), makeClone);
   funcFactory.addFunc("ConstantValue", new ConstantValue(), makeClone);
   funcFactory.addFunc("SpatialMap", new SpatialMap(), makeClone);
   
// Use either OptEM or classic Likelihoood.
   LogLike * logLike;
   bool useOptEM = params.bool_par("Use_OptEM");
   if (useOptEM) {
      logLike = new OptEM();
   } else {
      logLike = new LogLike();
   }

// Read in the Source model.
   std::string sourceModel = params.string_par("Source_model_file");
   logLike->readXml(sourceModel, funcFactory);
   
// Read in the Event data.
   std::string eventFile = params.string_par("event_file");
   int eventFileHdu = params.long_par("event_file_hdu");
   std::vector<std::string> eventFiles;
   resolve_fits_files(eventFile, eventFiles);
   std::vector<std::string>::const_iterator evIt = eventFiles.begin();
   for ( ; evIt != eventFiles.end(); evIt++) {
      logLike->getEvents(*evIt, eventFileHdu);
   }

// Compute the Event responses to the diffuse components.
   logLike->computeEventResponses();

// Set the verbosity level and convergence tolerance.
   int verbose = static_cast<int>(params.long_par("fit_verbosity"));
   double tol = params.double_par("fit_tolerance");
   std::vector<double> errors;

   if (useOptEM) {
      dynamic_cast<OptEM *>(logLike)->findMin(verbose);
   } else {
// Select an optimizer.
      optimizers::Optimizer *myOpt = 0;
      std::string optimizer = params.string_par("optimizer");
      if (optimizer == "LBFGS") {
         myOpt = new optimizers::Lbfgs(*logLike);
      } else if (optimizer == "MINUIT") {
         myOpt = new optimizers::Minuit(*logLike);
      } else if (optimizer == "DRMNGB") {
         myOpt = new optimizers::Drmngb(*logLike);
      }

// Do the fit.
      try {
         myOpt->find_min(verbose, tol);
      } catch (optimizers::Exception &eObj) {
         std::cerr << eObj.what() << std::endl;
      }

// Evaluate the uncertainties, if available.
      if (optimizer == "MINUIT") {
         errors = dynamic_cast<optimizers::Minuit *>(myOpt)->getUncertainty();
      } else if (optimizer == "DRMNGB") {
         int retCode = dynamic_cast<optimizers::Drmngb *>(myOpt)->getRetCode();
         std::cerr << "Drmngb return code: " << retCode;
         errors = dynamic_cast<optimizers::Drmngb *>(myOpt)->getUncertainty();
      } 
      delete myOpt;
   }
   print_fit_results(*logLike, errors);

// Write the model to the output xml file.
   std::string xmlFile = params.string_par("Source_model_output_file");
   std::string funcFileName = params.string_par("Function_models_file_name");

   std::cout << "Writing fitted model to " << xmlFile << std::endl;
   logLike->writeXml(xmlFile, funcFileName);

// Write the model to a flux-style output file.
   std::string xml_fluxFile = params.string_par("flux_style_model_file");
   if (xml_fluxFile != "none") {
      std::cout << "Writing flux-style xml model file to "
                << xml_fluxFile << std::endl;
      logLike->write_fluxXml(xml_fluxFile);
   }
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

void resolve_fits_files(std::string filename, 
                        std::vector<std::string> &files) {

   facilities::Util::expandEnvVar(&filename);
   files.clear();

   size_t pos = filename.find(".fits");
   if (pos != std::string::npos) {
// filename is the name of a FITS file. Return that as the sole
// element in the files vector.
      files.push_back(filename);
      return;
   } else {
// filename contains a list of fits files.
      readLines(filename, files);
      return;
   }
}

void readLines(std::string inputFile, std::vector<std::string> &lines) {

   facilities::Util::expandEnvVar(&inputFile);

   std::ifstream file(inputFile.c_str());
   lines.clear();
   std::string line;
   while (std::getline(file, line, '\n')) {
      lines.push_back(line);
   }
}
