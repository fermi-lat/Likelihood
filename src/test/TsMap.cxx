/**
 * @file TsMap.cxx
 * @brief Prototype standalone application for producing
 * "test-statistic" maps.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/test/TsMap.cxx,v 1.4 2003/11/12 22:01:40 jchiang Exp $
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
#include "optimizers/../src/PowerLaw.h"

#include "latResponse/IrfsFactory.h"

#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/SourceModel.h" 
#include "Likelihood/Source.h"
#include "Likelihood/ScData.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/SpatialMap.h"
#include "Likelihood/LogLike.h"
#include "Likelihood/RunParams.h"

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
   ExposureMap::readExposureFile(exposureFile);

// Create the response functions.
   std::string responseFuncs;
   params.getParam("Response_functions", responseFuncs);
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

// Add the prototypes for modeling spatial distributions.
   bool makeClone(false);
   funcFactory.addFunc("SkyDirFunction", new SkyDirFunction(), makeClone);
   funcFactory.addFunc("SpatialMap", new SpatialMap(), makeClone);
   
// Read in the Source model.
   std::string sourceModel;
   params.getParam("Source_model_file", sourceModel);
   logLike.readXml(sourceModel, funcFactory);
   
// Read in the Event data.
   std::string eventFile;
   params.getParam("event_file", eventFile);
   long eventFileHdu;
   params.getParam("event_file_hdu", eventFileHdu);
   std::vector<std::string> eventFiles;
   RunParams::resolve_fits_files(eventFile, eventFiles);
   std::vector<std::string>::const_iterator evIt = eventFiles.begin();
   for ( ; evIt != eventFiles.end(); evIt++) {
      logLike.getEvents(*evIt, eventFileHdu);
   }

// Compute the Event responses to the diffuse components.
   logLike.computeEventResponses();

// Select an optimizer.
   std::string optimizer;
   params.getParam("optimizer", optimizer);
   optimizers::Optimizer *myOpt = 0;
   if (optimizer == "LBFGS") {
      myOpt = new optimizers::Lbfgs(logLike);
   } else if (optimizer == "MINUIT") {
      myOpt = new optimizers::Minuit(logLike);
   } else if (optimizer == "DRMNGB") {
      myOpt = new optimizers::Drmngb(logLike);
   }

// Set the verbosity level and convergence tolerance.
   long verbose;
   params.getParam("fit_verbosity", verbose);
   double tol;
   params.getParam("fit_tolerance", tol);

// Retrieve map dimensions.
   std::string coordSystem;
   params.getParam("Coordinate_system", coordSystem);
   double lonMax;
   params.getParam("Longitude_max", lonMax);
   double lonMin;
   params.getParam("Longitude_min", lonMin);
   long nlon;
   params.getParam("Number_of_longitude_points", nlon);
   double latMax;
   params.getParam("Latitude_max", latMax);
   double latMin;
   params.getParam("Latitude_min", latMin);
   long nlat;
   params.getParam("Number_of_latitude_points", nlat);

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
   
   std::string outputFile;
   params.getParam("TS_map_file", outputFile);
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
   optimizers::PowerLaw *pl = new optimizers::PowerLaw(1, -2, 100.);
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
