/**
 * @file TsMap.cxx
 * @brief Prototype standalone application for producing
 * "test-statistic" maps.
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

#include "fitsio.h"

#include "optimizers/dArg.h"
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
void makeDoubleVector(double xmin, double xmax, int nx,
                      std::vector<double> &xVals);
void setPointSourceSpectrum(PointSource &src);
void write_fits_file(const std::string &filename, 
                     std::vector<double> &lon,
                     std::vector<double> &lat,
                     std::vector< std::vector<double> > &map,
                     const std::string &coordSystem);
void fitsReportError(FILE *stream, int status);

int main(int iargc, char* argv[]) {

// Read in the command-line parameters using HOOPS
   std::string filename("TsMap.par");
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
   if (responseFuncs == "COMBINED") {
      ResponseFunctions::addRespPtr(4, 
                                    irfsFactory.create("Glast25::Combined"));
   } else if (responseFuncs == "FRONT/BACK") {
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
   if (optimizer == "LBFGS") {
      myOpt = new optimizers::Lbfgs(logLike);
   } else if (optimizer == "MINUIT") {
      myOpt = new optimizers::Minuit(logLike);
   } else if (optimizer == "DRMNGB") {
      myOpt = new optimizers::Drmngb(logLike);
   }

// Set the verbosity level and convergence tolerance.
   int verbose = static_cast<int>(params.long_par("fit_verbosity"));
   double tol = params.double_par("fit_tolerance");

// Retrieve map dimensions.
   std::string coordSystem = params.string_par("Coordinate_system");
   double lonMax = params.double_par("Longitude_max");
   double lonMin = params.double_par("Longitude_min");
   int nlon = params.long_par("Number_of_longitude_points");
   double latMax = params.double_par("Latitude_max");
   double latMin = params.double_par("Latitude_min");
   int nlat = params.long_par("Number_of_latitude_points");

   std::vector<double> lonValues;
   makeDoubleVector(lonMin, lonMax, nlon, lonValues);
   std::vector<double> latValues;
   makeDoubleVector(latMin, latMax, nlat, latValues);

   std::string srcName("testSource");
   optimizers::dArg dummy(1.);

// Save the best-fit value in the null hypothesis.
   myOpt->find_min(verbose, tol);
   double logLike0 = logLike(dummy);

// Create a putative PointSource to be moved to each map location.
   PointSource testSrc;
   setPointSourceSpectrum(testSrc);
   testSrc.setName(srcName);

   bool computeExposure(true);

   std::vector< std::vector<double> > myMap(lonValues.size());

// Loop over map positions.
   for (unsigned int ii = 0; ii < lonValues.size(); ii++) {
      myMap[ii].reserve(latValues.size());
      for (unsigned int jj = 0; jj < latValues.size(); jj++) {
         if (coordSystem == "CEL") {
            testSrc.setDir(lonValues[ii], latValues[jj], computeExposure);
         } else {
            testSrc.setGalDir(lonValues[ii], latValues[jj], computeExposure);
         }
         logLike.addSource(&testSrc);
         myOpt->find_min(verbose, tol);
         myMap[ii].push_back(2.*(logLike(dummy) - logLike0));
         if (verbose > 0) {
            std::cout << lonValues[ii] << "  "
                      << latValues[jj] << "  "
                      << myMap[ii][jj] 
                      << std::endl;
         }
         logLike.deleteSource(srcName);
      }
   }
   
   std::string outputFile = params.string_par("TS_map_file");
   write_fits_file(outputFile, lonValues, latValues, myMap, coordSystem);

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

void makeDoubleVector(double xmin, double xmax, int nx,
                      std::vector<double> &xVals) {
   xVals.clear();
   xVals.reserve(nx);
   double xstep = (xmax - xmin)/(nx-1);
   for (int i = 0; i < nx; i++) {
      xVals.push_back(xstep*i + xmin);
   }
}

void setPointSourceSpectrum(PointSource &src) {
   PowerLaw *pl = new PowerLaw(1, -2, 100.);
   optimizers::Parameter indexParam = pl->getParam("Index");
   indexParam.setBounds(-3.5, -1.);
   pl->setParam(indexParam);
   optimizers::Parameter prefactorParam = pl->getParam("Prefactor");
   prefactorParam.setBounds(1e-10, 1e3);
   prefactorParam.setScale(1e-9);
   pl->setParam(prefactorParam);
   src.setSpectrum(pl);
}

void write_fits_file(const std::string &filename, 
                     std::vector<double> &lon,
                     std::vector<double> &lat,
                     std::vector< std::vector<double> > &map,
                     const std::string &coordSystem) {

   fitsfile *fptr;
   int status = 0;
   
// Always overwrite an existing file.
   remove(filename.c_str());
   fits_create_file(&fptr, filename.c_str(), &status);
   fitsReportError(stderr, status);

   int bitpix = DOUBLE_IMG;
   long naxis = 2;
   long naxes[] = {lon.size(), lat.size()};
   fits_create_img(fptr, bitpix, naxis, naxes, &status);
   fitsReportError(stderr, status);

// Repack exposure into a C array.
   double *map_array = new double[lon.size()*lat.size()];
   int indx = 0;
   for (unsigned int i = 0; i < lon.size(); i++) {
      for (unsigned int j = 0; j < lat.size(); j++) {
         map_array[indx] = map[i][j];
         indx++;
      }
   }

// Write the exposure map data.
   long group = 0;
   long firstelem = 1;
   long nelements = lon.size()*lat.size();
   fits_write_img_dbl(fptr, group, firstelem, nelements,
                      map_array, &status);
   delete[] map_array;
   fitsReportError(stderr, status);

// Write some keywords.
   double l0 = lon[0];
   fits_update_key(fptr, TDOUBLE, "CRVAL1", &l0, 
                   "longitude of reference pixel", &status);
   fitsReportError(stderr, status);
   double b0 = lat[0];
   fits_update_key(fptr, TDOUBLE, "CRVAL2", &b0, 
                   "latitude of reference pixel", &status);
   fitsReportError(stderr, status);
   
   double lstep = lon[1] - lon[0];
   fits_update_key(fptr, TDOUBLE, "CDELT1", &lstep, 
                   "longitude step at ref. pixel", &status);
   fitsReportError(stderr, status);
   double bstep = lat[1] - lat[0];
   fits_update_key(fptr, TDOUBLE, "CDELT2", &bstep, 
                   "latitude step at ref. pixel", &status);
   fitsReportError(stderr, status);
   
   float crpix1 = lstep/2.;
   fits_update_key(fptr, TFLOAT, "CRPIX1", &crpix1, 
                   "reference pixel for longitude coordinate", &status);
   fitsReportError(stderr, status);
   float crpix2 = bstep/2.;
   fits_update_key(fptr, TFLOAT, "CRPIX2", &crpix2, 
                   "reference pixel for latitude coordinate", &status);
   fitsReportError(stderr, status);
   
   if (coordSystem == "CEL") {
      char *ctype1 = "RA---CAR";
      fits_update_key(fptr, TSTRING, "CTYPE1", ctype1, 
                      "right ascension", &status);
      fitsReportError(stderr, status);
      char *ctype2 = "DEC--CAR";
      fits_update_key(fptr, TSTRING, "CTYPE2", ctype2, 
                      "declination", &status);
      fitsReportError(stderr, status);
   } else if (coordSystem == "GAL") {
      char *ctype1 = "GLON-CAR";
      fits_update_key(fptr, TSTRING, "CTYPE1", ctype1, 
                      "Galactic longitude", &status);
      fitsReportError(stderr, status);
      char *ctype2 = "GLAT-CAR";
      fits_update_key(fptr, TSTRING, "CTYPE2", ctype2, 
                      "Galactic latitude", &status);
      fitsReportError(stderr, status);
   }
   
   fits_close_file(fptr, &status);
   fitsReportError(stderr, status);
   
   return;
}

void fitsReportError(FILE *stream, int status) {
   fits_report_error(stream, status);
   if (status != 0) {
      throw std::string("writeExposureFile: cfitsio error.");
   }
}
