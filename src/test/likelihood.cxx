/**
 * @file likelihood.cxx
 * @brief Prototype standalone application for the Likelihood tool.
 * @author J. Chiang
 *
 * $Header$
 */

#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>

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
#include "Likelihood/RunParams.h"
#include "PowerLaw.h"
#include "Gaussian.h"
#include "AbsEdge.h"

using namespace Likelihood;

void print_fit_results(SourceModel &stat);

int main(int iargc, char* argv[]) {

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
   ScData::readData(scFile, scHdu);

// Read in the exposure map file.
   std::string exposureFile = params.string_par("Exposure_map_file");
   ExposureMap::readExposureFile(exposureFile);

// Create the response functions.
   std::string responseFuncs = params.string_par("Response_functions");
   latResponse::IrfsFactory irfsFactory;
   if (responseFuncs == std::string("Combined")) {
      ResponseFunctions::addRespPtr(4, 
                                    irfsFactory.create("Glast25::Combined"));
   } else if (responseFuncs == std::string("Front/Back")) {
      ResponseFunctions::addRespPtr(2, irfsFactory.create("Glast25::Front"));
      ResponseFunctions::addRespPtr(3, irfsFactory.create("Glast25::Back"));
   }

// Use unbinned log-likelihood as the objective function.
   LogLike logLike;

// Fill a FunctionFactory with Function object prototypes for source
// modeling.
   optimizers::FunctionFactory funcFactory;

// Add the standard prototypes for modeling spectra,
   bool makeClone(false);
   funcFactory.addFunc("PowerLaw", new PowerLaw(), makeClone);
   funcFactory.addFunc("Gaussian", new Gaussian(), makeClone);
   funcFactory.addFunc("AbsEdge", new AbsEdge(), makeClone);

// and the prototypes for modeling spatial distributions.
   funcFactory.addFunc("SkyDirFunction", new SkyDirFunction(), makeClone);
   funcFactory.addFunc("ConstantValue", new ConstantValue(), makeClone);
   funcFactory.addFunc("SpatialMap", new SpatialMap(), makeClone);
   
// Read in the Source model.
   std::string sourceModel = params.string_par("Source_model_file");
   logLike.readXml(sourceModel, funcFactory);
   
// Read in the Event data.
   std::string eventFile = params.string_par("event_file");
   int eventFileHdu = params.long_par("event_file_hdu");
   logLike.getEvents(eventFile, eventFileHdu);

// Compute the Event responses to the diffuse components.
   logLike.computeEventResponses();

// Select an optimizer.
   std::string optimizer = params.string_par("optimizer");
   optimizers::Optimizer *myOpt;
   if (optimizer == "Lbfgs") {
      myOpt = new optimizers::Lbfgs(logLike);
   } else if (optimizer == "Minuit") {
      myOpt = new optimizers::Minuit(logLike);
   } else if (optimizer == "Drmngb") {
      myOpt = new optimizers::Drmngb(logLike);
   }

// Set the verbosity level and convergence tolerance.
   int verbose = static_cast<int>(params.long_par("fit_verbosity"));
   double tol = params.double_par("fit_tolerance");

// Do the fit.
   myOpt->find_min(verbose, tol);

   print_fit_results(logLike);

// Write the model to the output xml file.
   std::string xmlFile = params.string_par("Source_model_output_file");
   std::string funcFileName = params.string_par("Function_models_file_name");

   std::cout << "Writing fitted model to " << xmlFile << std::endl;
   logLike.writeXml(xmlFile, funcFileName);

   delete myOpt;
}

void print_fit_results(SourceModel &stat) {
   std::vector<std::string> srcNames;
   stat.getSrcNames(srcNames);
   std::vector<optimizers::Parameter> parameters;
   for (unsigned int i = 0; i < srcNames.size(); i++) {
      Source *src = stat.getSource(srcNames[i]);
      Source::FuncMap srcFuncs = src->getSrcFuncs();
      srcFuncs["Spectrum"]->getParams(parameters);
      std::cout << "\n" << srcNames[i] << ":\n";
      for (unsigned int i = 0; i < parameters.size(); i++)
         std::cout << parameters[i].getName() << ": "
                   << parameters[i].getValue() << std::endl;
      std::cout << "Npred: "
                << src->Npred() << std::endl;
   }
}
