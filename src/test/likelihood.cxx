#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>

#include "astro/SkyDir.h"

#include "latResponse/../src/Table.h"
#include "latResponse/IrfsFactory.h"

#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/SourceModel.h" 
#include "Likelihood/Event.h"
#include "Likelihood/Source.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/ScData.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/SpectrumFactory.h"
#include "Likelihood/SourceFactory.h"
#include "Likelihood/FitsImage.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/SpatialMap.h"
#include "Likelihood/ConstantValue.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/LogLike.h"
#include "Likelihood/RunParams.h"
#include "PowerLaw.h"
#include "Gaussian.h"
#include "AbsEdge.h"

#include "optimizers/Function.h"
#include "optimizers/SumFunction.h"
#include "optimizers/ProductFunction.h"
#include "optimizers/FunctionFactory.h"
#include "optimizers/Lbfgs.h"
#include "optimizers/Minuit.h"
#include "optimizers/Exception.h"

using namespace Likelihood;

void print_fit_results(SourceModel &stat);

int main(int iargc, char* argv[]) {

   std::string filename("likelihood.par");
   delete argv[0];
   argv[0] = strdup(filename.c_str());

   RunParams params(iargc, argv);

   std::string roiCutsFile = params.string_par("ROI_cuts_file");
   RoiCuts::setCuts(roiCutsFile);

   std::string scFile = params.string_par("Spacecraft_file");
   int scHdu = static_cast<int>(params.long_par("Spacecraft_file_hdu"));
   ScData::readData(scFile, scHdu);

   std::string exposureFile = params.string_par("Exposure_map_file");
   ExposureMap::readExposureFile(exposureFile);

   std::string responseFuncs = params.string_par("Response_functions");
   latResponse::IrfsFactory irfsFactory;
   if (responseFuncs == std::string("Combined")) {
      ResponseFunctions::addRespPtr(4, 
                                    irfsFactory.create("Glast25::Combined"));
   } else if (responseFuncs == std::string("Front/Back")) {
      ResponseFunctions::addRespPtr(2, irfsFactory.create("Glast25::Front"));
      ResponseFunctions::addRespPtr(3, irfsFactory.create("Glast25::Back"));
   }

   LogLike logLike;

   std::string sourceModel = params.string_par("Source_model_file");
   optimizers::FunctionFactory funcFactory;

// Add standard prototypes for modeling spectra,
   bool makeClone(false);
   funcFactory.addFunc("PowerLaw", new PowerLaw(), makeClone);
   funcFactory.addFunc("Gaussian", new Gaussian(), makeClone);
   funcFactory.addFunc("AbsEdge", new AbsEdge(), makeClone);

// and some prototypes for modeling spatial distributions.
   funcFactory.addFunc("SkyDirFunction", new SkyDirFunction(), makeClone);
   funcFactory.addFunc("ConstantValue", new ConstantValue(), makeClone);
   funcFactory.addFunc("SpatialMap", new SpatialMap(), makeClone);
   
   logLike.readXml(sourceModel, funcFactory);
   
   std::string eventFile = params.string_par("event_file");
   int eventFileHdu = params.long_par("event_file_hdu");
   logLike.getEvents(eventFile, eventFileHdu);

   logLike.computeEventResponses();

   std::string optimizer = params.string_par("optimizer");
   optimizers::Optimizer *myOpt;
   if (optimizer == "Lbfgs") {
      myOpt = new optimizers::Lbfgs(logLike);
   } else if (optimizer == "Minuit") {
      myOpt = new optimizers::Minuit(logLike);
   }
//  else if (optimizer == "Drmngb") {
//       myOpt = new optimizers::Drmngb(logLike);
//    }

   int verbose = static_cast<int>(params.long_par("verbosity"));
   myOpt->find_min(verbose, 1e-3);

   print_fit_results(logLike);

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
