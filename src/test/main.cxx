// test program for Likelihood

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include "astro/SkyDir.h"

//  include everything for the compiler to test

#include "Likelihood/SourceModel.h" 
#include "Likelihood/Table.h"
#include "Likelihood/Event.h"
#include "Likelihood/Source.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/Response.h"
#include "Likelihood/Aeff.h"
#include "Likelihood/Psf.h"
#include "Likelihood/ScData.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/SpectrumFactory.h"
#include "Likelihood/SourceFactory.h"
#include "Likelihood/FitsImage.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/SpatialMap.h"
#include "Likelihood/ConstantValue.h"
#include "Likelihood/DiffuseSource.h"
#include "logLike_ptsrc.h"
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

// Various using declarations/directives for testing purposes only
using namespace Likelihood;
using optimizers::Parameter;
using optimizers::Function;
using optimizers::SumFunction;
using optimizers::ProductFunction;
using optimizers::FunctionFactory;
using optimizers::dArg;
using optimizers::Lbfgs;
using optimizers::Minuit;

void read_SC_Response_data();
void test_SourceModel_class();
void report_SrcModel_values(const SourceModel &SrcModel);
void test_Table_class();
void test_Event_class();
void test_PointSource_class();
void test_Aeff_class();
void test_Psf_class();
void test_logLike_ptsrc();
void fit_3C279();
void fit_anti_center();
void test_SpectrumFactory();
void test_SourceFactory();
void test_FitsImage();
void test_ExposureMap();
void test_SpatialMap();
void test_DiffuseSource();
void fit_DiffuseSource();
void print_fit_results(SourceModel &stat);
void test_cfitsio();
void test_FunctionFactory();

std::string root_path;
std::string test_path;

int main() {
   read_SC_Response_data();
//    test_SourceModel_class();
//    test_Table_class();
//    test_Event_class();
//    test_PointSource_class();
//    test_Aeff_class();
//    test_Psf_class();
//    test_logLike_ptsrc();
//    test_SpectrumFactory();
//    fit_3C279();
//    fit_anti_center();
//    test_FitsImage();
//    test_ExposureMap();
//    test_SpatialMap();
//    test_DiffuseSource();
//    test_FunctionFactory();
//    test_SourceFactory();
   fit_DiffuseSource();
   return 0;
}

void test_FunctionFactory() {
   std::cout << "*** test FunctionFactory ***" << std::endl;

   FunctionFactory funcFactory;

   try {
// Add in the Functions for modeling spectra...
      funcFactory.addFunc("PowerLaw", new PowerLaw(), false);
      funcFactory.addFunc("Gaussian", new Gaussian(), false);
      funcFactory.addFunc("AbsEdge", new AbsEdge(), false);

// and the Functions for spatial modeling.
      funcFactory.addFunc("SkyDirFunction", new SkyDirFunction(), false);
      funcFactory.addFunc("ConstantValue", new ConstantValue(), false);
      funcFactory.addFunc("SpatialMap", new SpatialMap(), false);

   } catch (optimizers::Exception &eObj) {
      std::cout << eObj.what() << std::endl;
   }

// Read in the prototypes from the XML file.
   std::string xmlFile = root_path + "/xml/A1_Functions.xml";

   try {
      funcFactory.readXml(xmlFile);
   } catch (optimizers::Exception &eObj) {
      std::cout << eObj.what() << std::endl;
   }

   std::cout << "*** FunctionFactory tests completed ***" << std::endl;
}

void fit_DiffuseSource() {
   std::cout << "*** fit_DiffuseSource ***" << std::endl;

// // center the ROI on 3C 279
//    double ra0 = 193.98;
//    double dec0 = -5.82;

//    RoiCuts::setCuts(ra0, dec0, 20.);

   std::string roiFile = root_path + "/xml/RoiCuts.xml";
   RoiCuts::setCuts(roiFile);

   double ra0, dec0;
   RoiCuts::getRaDec(ra0, dec0);
   std::cout << "ROI center: " << ra0 << "  " << dec0;

// root name for the observation data files
   std::string obs_root = "diffuse_test_5";

// read in the spacecraft data
   std::string sc_file = test_path + "Data/" + obs_root + "_sc_0000";
   int sc_hdu = 2;
   ScData::readData(sc_file, sc_hdu);

   std::string expfile = test_path + "Data/exp_" + obs_root + "_new.fits";
//   std::string expfile = "exp_" + obs_root + ".fits";

// compute a new exposure map for these data
//   ExposureMap::computeMap(expfile, 30., 60, 60, 10);

// must read in the exposure file prior to creating the SourceFactory
// object since it contains DiffuseSources
   ExposureMap::readExposureFile(expfile);

// Create the FunctionFactory and SourceFactory.
   optimizers::FunctionFactory funcFactory;
   try {
// Add in the Functions for modeling spectra...
      funcFactory.addFunc("PowerLaw", new PowerLaw(), false);
      funcFactory.addFunc("Gaussian", new Gaussian(), false);
      funcFactory.addFunc("AbsEdge", new AbsEdge(), false);

// and the Functions for spatial modeling.
      funcFactory.addFunc("SkyDirFunction", new SkyDirFunction(), false);
      funcFactory.addFunc("ConstantValue", new ConstantValue(), false);
      funcFactory.addFunc("SpatialMap", new SpatialMap(), false);

   } catch (optimizers::Exception &eObj) {
      std::cout << eObj.what() << std::endl;
   }

   SourceFactory srcFactory;

// Read in the prototypes from the XML file.
   std::string xmlFile = root_path + "/xml/A1_Sources.xml";

   srcFactory.readXml(xmlFile, funcFactory);

   DiffuseSource *ourGalaxy = dynamic_cast<DiffuseSource *>
      (srcFactory.create("Galactic Diffuse Emission"));

   std::cout << ourGalaxy->getName() << std::endl;

   DiffuseSource *extragalactic = dynamic_cast<DiffuseSource *>
      (srcFactory.create("Extragalactic Diffuse Emission"));

   std::cout << extragalactic->getName() << std::endl;

// 3C 279
   Source *_3c279 = srcFactory.create("Bright Point Source");
   _3c279->setDir(ra0, dec0);
   _3c279->setName("3C 279");

   std::cout << _3c279->getName() << std::endl;

// create the Statistic
   logLike_ptsrc logLike;

// add the Sources
   logLike.addSource(ourGalaxy);
   logLike.addSource(extragalactic);
   logLike.addSource(_3c279);

// read in the data
   std::string event_file = test_path + "Data/" + obs_root + "_0000";
   logLike.getEvents(event_file, 2);

// There are a few options for computing the DiffuseSource Event responses:

// individually...
//     logLike.computeEventResponses(*ourGalaxy);
//     logLike.computeEventResponses(*extragalactic);

// by constructing a vector of the targeted DiffuseSources...
//     std::vector<DiffuseSource> srcs;
//     srcs.push_back(*ourGalaxy);
//     srcs.push_back(*extragalactic);
//     logLike.computeEventResponses(srcs);

// or the default way, for all of the DiffuseSources in the SourceModel...
   logLike.computeEventResponses();

// Derivative tests:
   dArg dummy(1.);
   std::vector<double> derivs;
   std::vector<double> params;

   double statval0 = logLike(dummy);

   try {
      logLike.getFreeDerivs(dummy, derivs);
      logLike.getFreeParamValues(params);
   } catch (Exception &eObj) {
      std::cout << eObj.what() << std::endl;
   } catch (optimizers::Exception &eObj) {
      std::cout << eObj.what() << std::endl;
   }
      
   try {
      for (unsigned int i = 0; i < params.size(); i++) {
         std::vector<double> new_params = params;
         double dparam = params[i]*1e-7;
         new_params[i] += dparam;
         logLike.setFreeParamValues(new_params);
         double statval = logLike(dummy);
         std::cout << derivs[i] << "  " 
                   << (statval - statval0)/dparam << std::endl;
      }
   } catch (Exception &eObj) {
      std::cout << eObj.what() << std::endl;
   } catch (optimizers::Exception &eObj) {
      std::cout << eObj.what() << std::endl;
   }

   std::vector<std::string> srcNames;
   logLike.getSrcNames(srcNames);
   
   int verbose = 3;

   Lbfgs myOpt(logLike);
   try {
      myOpt.find_min(verbose);
   } catch (Exception eObj) {
      std::cerr << eObj.what() << std::endl;
   } catch (optimizers::Exception eObj) {
      std::cerr << eObj.what() << std::endl;
   }
   print_fit_results(logLike);
      
// replace (or add) each Source in srcFactory for later use
   std::vector<std::string> factoryNames;
   srcFactory.fetchSrcNames(factoryNames);
   for (unsigned int i = 0; i < srcNames.size(); i++) {
      std::cerr << "Adding " << srcNames[i] 
                << " to srcFactory." << std::endl;
      Source *src = logLike.getSource(srcNames[i]);
      srcFactory.replaceSource(src);
   }

// try to fit again using srcFactory Sources and the MINUIT optimizer
   logLike.deleteAllSources();

   for (unsigned int i = 0; i < srcNames.size(); i++) {
      Source *src = srcFactory.create(srcNames[i]);
      logLike.addSource(src);
   }

   verbose = 3;
   Minuit myMinuitObj(logLike);
   try {
      myMinuitObj.find_min(verbose, .0001);
   } catch (Exception eObj) {
      std::cerr << eObj.what() << std::endl;
   } catch (optimizers::Exception eObj) {
      std::cerr << eObj.what() << std::endl;
   }

   std::vector<double> sig = myMinuitObj.getUncertainty();
   for (unsigned int i=0; i < sig.size(); i++) {
      std::cout << i << "  " << sig[i] << std::endl;
   }
   print_fit_results(logLike);

// Write out the fitted model as an xml file.
   logLike.writeXml("fittedModel.xml");

   std::cout << "*** fit_DiffuseSource: all tests completed ***\n" 
             << std::endl;

} // fit_DiffuseSource

void test_DiffuseSource() {

   std::cout << "*** test_DiffuseSource ***" << std::endl;

// center the ROI on the Galactic Center
//     double ra0 = 266.404;
//     double dec0 = -28.9364;
   double ra0 = 193.98;
   double dec0 = -5.82;
   RoiCuts::setCuts(ra0, dec0, 20.);

//  /* read in the spacecraft data */
//     std::string sc_file = test_path + "Data/one_src_sc_0000";
//     int sc_hdu = 2;
//     ScData::readData(sc_file, sc_hdu);

   std::string expfile = test_path + "Data/exp_diffuse_test_5_new.fits";
   ExposureMap::readExposureFile(expfile);

   std::string galfile = test_path + "Data/gas.cel";
   SpatialMap galacticModel(galfile);

   DiffuseSource ourGalaxy(&galacticModel);
   ourGalaxy.setName("Milky Way");

// Provide ourGalaxy with a power-law spectrum.
// Since gas.cel is integral photon flux above 100 MeV, the power-law
// Prefactor needs also to be scaled by the differential flux at 1 MeV.
   PowerLaw gal_pl(pow(100., -2.1)*1.1*pow(100., 1.1), -2.1, 100.);
   ourGalaxy.setSpectrum(&gal_pl);

// Output the fluxDensity along a meridian passing through the
// position of 3C 279.  Unfortunately, one must do this by creating
// Event objects, which entail executing the computeResponse() method
// for each Event. One must also assume an attitude for the LAT.  Here
// we simple take the z-axis to be in the direction of 3C 279.

   std::vector<Event> my_Events;
   int nevents = 20;
   double ra = ra0;
   double decmin = dec0 - 5.;
   double decstep = 2.*(dec0 - decmin)/(nevents-1.);
   double energy = 100;
   double time = 60;
   double muZenith = -1; // this is presently irrelevant but must be provided
   for (int i = 0; i < nevents; i++) {
      double dec = decstep*i + decmin;
      my_Events.push_back(Event(ra, dec, energy, time, ra, dec, muZenith));
      my_Events[i].computeResponse(ourGalaxy);

      std::cout << dec << "  "
                << ourGalaxy.fluxDensity(my_Events[i]) << "  "
                << ourGalaxy.Npred()
                << std::endl;
   }

   std::cout << "*** test_DiffuseSource: all tests completed ***\n" 
             << std::endl;

} // test_DiffuseSource

void test_SpatialMap() {

   std::cout << "*** test_SpatialMap ***" << std::endl;

   std::string fitsfile = test_path + "Data/gas.cel";
   SpatialMap my_map(fitsfile);

   double dec = 0;

   std::vector<double> ra;
   for (double ra_val = -180; ra_val < 180.; ra_val += 5)
      ra.push_back(ra_val);

   for (unsigned int i = 0; i < ra.size(); i++) {
      astro::SkyDir dir(ra[i], dec);
      SkyDirArg sArg(dir);
      std::cout << ra[i] << "  "
                << my_map(sArg) << std::endl;
   }

   std::cout << "*** test_SpatialMap: all tests completed ***\n" << std::endl;

} // test_SpatialMap

void test_ExposureMap() {
   std::cout << "*** test_ExposureMap ***" << std::endl;

   std::string fitsfile = test_path + "Data/test.fits";
   ExposureMap::readExposureFile(fitsfile);

   ExposureMap *emap = ExposureMap::instance();

   std::vector<double> energies;
   emap->fetchEnergies(energies);
   for (unsigned int i = 0; i < energies.size(); i += energies.size()/5)
      std::cout << energies[i] << " ";
   std::cout << std::endl;

   std::vector< std::valarray<double> > exposure;
   emap->fetchExposure(exposure);
   for (unsigned int i = 0; i < exposure.size(); i += exposure.size()/5) {
      for (unsigned int j = 0; j < exposure[i].size(); 
           j += exposure[i].size()/10) {
         std::cout << exposure[i][j] << " ";
      }
      std::cout << std::endl;
   }

   std::cout << "*** test_ExposureMap: all tests completed ***\n" << std::endl;

} // test_ExposureMap

void test_FitsImage() {

   std::cout << "*** test_FitsImage ***" << std::endl;

//   std::string fitsfile = test_path + "Data/anti_center_exposure.fits";
   std::string fitsfile = test_path + "Data/test.fits";
   try {
      FitsImage expImage(fitsfile);

   std::cerr << "AxisDims: ";
   std::vector<int> axisDims;
   expImage.fetchAxisDims(axisDims);
   for (unsigned int i = 0; i < axisDims.size(); i++) 
      std::cerr << axisDims[i] << "  ";
   std::cerr << std::endl;
   std::cerr << std::endl;

   std::cerr << "AxisNames: ";
   std::vector<std::string> axisNames;
   expImage.fetchAxisNames(axisNames);
   for (unsigned int i = 0; i < axisNames.size(); i++)
      std::cout << axisNames[i] << "  ";
   std::cout << std::endl;
   std::cout << std::endl;

   std::cerr << "AxisVectors: \n";
   std::vector< std::vector<double> > axisVectors;
   axisVectors.resize(axisNames.size());
   for (unsigned int i = 0; i < axisNames.size(); i++)
      expImage.fetchAxisVector(i, axisVectors[i]);
   for (unsigned int i = 0; i < axisVectors.size(); i++) {
      for (unsigned int j = 0; j < axisVectors[i].size(); 
           j += axisVectors[i].size()/5)
         std::cout << axisVectors[i][j] << "  ";
      std::cout << std::endl;
   }
   std::cout << std::endl;

   std::cerr << "CelestialArrays: \n";
   std::valarray<double> lonArray;
   std::valarray<double> latArray;
   expImage.fetchCelestialArrays(lonArray, latArray);
   for (unsigned int i = 0; i < lonArray.size(); i += lonArray.size()/7)
      std::cout << lonArray[i] << "  ";
   std::cout << std::endl;
   for (unsigned int i = 0; i < latArray.size(); i += latArray.size()/7)
      std::cout << latArray[i] << "  ";
   std::cout << std::endl;
   
   std::valarray<double> imageData;
   expImage.fetchImageData(imageData);
   } catch (Exception eObj) {
      std::cerr << eObj.what() << std::endl;
   } 
   catch (...) {
   }

   std::cout << "*** test_FitsImage: all tests completed ***\n" << std::endl;
   
} // test_FitsImage

void test_SourceFactory() {

   std::cout << "*** test_SourceFactory ***" << std::endl;

// center the ROI on 3C 279
   double ra0 = 193.98;
   double dec0 = -5.82;

   RoiCuts::setCuts(ra0, dec0, 20.);

// root name for the observation data files
   std::string obs_root = "diffuse_test_5";

// read in the spacecraft data
   std::string sc_file = test_path + "Data/" + obs_root + "_sc_0000";
   int sc_hdu = 2;
   ScData::readData(sc_file, sc_hdu);

   std::string expfile = test_path + "Data/exp_" + obs_root + "_new.fits";
//   std::string expfile = "exp_" + obs_root + ".fits";

// compute a new exposure map for these data
//   ExposureMap::computeMap(expfile, 30., 60, 60, 10);

// must read in the exposure file prior to creating the SourceFactory
// object since it contains DiffuseSources
   ExposureMap::readExposureFile(expfile);

// Create the FunctionFactory and SourceFactory.
   optimizers::FunctionFactory funcFactory;
   try {
// Add in the Functions for modeling spectra...
      funcFactory.addFunc("PowerLaw", new PowerLaw(), false);
      funcFactory.addFunc("Gaussian", new Gaussian(), false);
      funcFactory.addFunc("AbsEdge", new AbsEdge(), false);

// and the Functions for spatial modeling.
      funcFactory.addFunc("SkyDirFunction", new SkyDirFunction(), false);
      funcFactory.addFunc("ConstantValue", new ConstantValue(), false);
      funcFactory.addFunc("SpatialMap", new SpatialMap(), false);

   } catch (optimizers::Exception &eObj) {
      std::cout << eObj.what() << std::endl;
   }

   SourceFactory srcFactory;

// Read in the prototypes from the XML file.
   std::string xmlFile = root_path + "/xml/A1_Sources.xml";

   srcFactory.readXml(xmlFile, funcFactory);

   std::cout << "*** test_SourceFactory: all tests completed ***\n" 
             << std::endl;

} // test_SourceFactory

void test_SpectrumFactory() {

   std::cout << "*** test_SpectrumFactory ***" << std::endl;

   SpectrumFactory myFactory;

   Function *pl_continuum = myFactory.makeFunction("PowerLaw");
   double parray[] = {1., -2, 1.};
   std::vector<double> params(parray, parray+3);
   pl_continuum->setParamValues(params);

   Function *emission_line = myFactory.makeFunction("Gaussian");
   params[0] = 10.; params[1] = 1.; params[2] = 0.1;
   emission_line->setParamValues(params);

   SumFunction spectrum(*pl_continuum, *emission_line);

   Function *edge = myFactory.makeFunction("AbsEdge");
   params[0] = 5; params[1] = 1.2; params[2] = -3;
   edge->setParamValues(params);

   ProductFunction absorbed_spec(*edge, spectrum);

// add this new Function to myFactory

   myFactory.addFunc("absorbed_spec", &absorbed_spec);

// list the available guys

   myFactory.listFunctions();

// clone this CompositeFunction using myFactory and see if it works 
// as expected after tweaking a parameter

   Function *abs_spec = myFactory.makeFunction("absorbed_spec");

   abs_spec->getParamValues(params);
   params[0] *= 1.5;  //modify the AbsEdge depth Tau0
   abs_spec->setParamValues(params);

   int nmax = 100;
   double xmin = 0.1;
   double xmax = 10.;
   double xstep = log(xmax/xmin)/(nmax-1);
   for (int i = 0; i < nmax; i++) {
      double x = xmin*exp(xstep*i);
      dArg xarg(x);
      std::cout << x << "  " 
                << absorbed_spec(xarg) << "  "
                << (*abs_spec)(xarg) 
                << std::endl;
   }
   
// check the derivatives of absorbed_spec
   absorbed_spec.getParamValues(params);
   dArg xarg(1.5);
   double spec_value = absorbed_spec(xarg);
   std::vector<double> derivs;
   absorbed_spec.getDerivs(xarg, derivs);

   double eps = 1e-7;
   for (unsigned int i = 0; i < params.size(); i++) {
      std::vector<double> new_params = params;
      double dparam = eps*new_params[i];
      new_params[i] += dparam;
      absorbed_spec.setParamValues(new_params);
      double new_spec_value = absorbed_spec(xarg);
      double num_deriv = (new_spec_value - spec_value)/dparam;
      std::cerr << derivs[i] << "  "
                << num_deriv << "  "
                << derivs[i]/num_deriv << "  "
                << std::endl;
   }
   std::cout << "*** test_SpectrumFactory: all tests completed ***\n" 
             << std::endl;

}  // test_SpectrumFactory

void fit_anti_center() {

   std::cout << "*** fit_anti_center ***" << std::endl;

// set the ROI cuts
   RoiCuts::setCuts(83.57, 22.01, 20.);

/* read in the spacecraft data */
   std::string sc_file = test_path + "Data/anti_center_sc_0000";
   int sc_hdu = 2;
   ScData::readData(sc_file, sc_hdu);

// add the sources to the model

   logLike_ptsrc logLike;

// Crab Pulsar
   double ra = 83.57;
   double dec = 22.01;
   PointSource Crab(ra, dec);

//   double Prefactor = 6.46e-4;
//   double Prefactor = 6.46e1;
   double Prefactor = 27.;
   double Gamma = -2.19;
   double Escale = 100.;

   PowerLaw Crab_pl(Prefactor, Gamma, Escale);

//set limits on index
   Parameter indexParam = Crab_pl.getParam("Index");
   indexParam.setBounds(-3.5, -1.);
   Crab_pl.setParam(indexParam);
   
// set limits on normalization
   Parameter prefactorParam = Crab_pl.getParam("Prefactor");
   prefactorParam.setBounds(1e-3, 1e3);
   prefactorParam.setScale(1e-9);
   Crab_pl.setParam(prefactorParam);

   Crab.setSpectrum(&Crab_pl);
   Crab.setName("Crab Pulsar");
   logLike.addSource(&Crab);

// Geminga Pulsar
   ra = 98.49;
   dec = 17.86;
   PointSource Geminga(ra, dec);

//   Prefactor = 4.866e-5;
   Prefactor = 23.29;
   Gamma = -1.66;
   Escale = 100.;

   PowerLaw Geminga_pl(Prefactor, Gamma, Escale);
//set limits on index
   indexParam = Geminga_pl.getParam("Index");
   indexParam.setBounds(-3.5, -1.);
   Geminga_pl.setParam(indexParam);
// set limits on normalization
   prefactorParam = Geminga_pl.getParam("Prefactor");
   prefactorParam.setBounds(1e-3, 1e3);  
   prefactorParam.setScale(1e-9);
   Geminga_pl.setParam(prefactorParam);

   Geminga.setSpectrum(&Geminga_pl);
   Geminga.setName("Geminga Pulsar");
   logLike.addSource(&Geminga);

// PKS 0528+134
   ra = 82.74;
   dec = 13.38;
   PointSource _0528(ra, dec);

//   Prefactor = 1.135e-3;
   Prefactor = 13.6;
   Gamma = -2.46;
   Escale = 100.;

   PowerLaw _0528_pl(Prefactor, Gamma, Escale);
//set limits on index
   indexParam = _0528_pl.getParam("Index");
   indexParam.setBounds(-3.5, -1.);  
   _0528_pl.setParam(indexParam);
// set limits on normalization
   prefactorParam = _0528_pl.getParam("Prefactor");
   prefactorParam.setBounds(1e-3, 1e3);
   prefactorParam.setScale(1e-9);
   _0528_pl.setParam(prefactorParam);

   _0528.setSpectrum(&_0528_pl);
   _0528.setName("PKS 0528+134");
   logLike.addSource(&_0528);

// read in the event data
   std::string event_file = test_path + "Data/anti_center_0000";
   logLike.getEvents(event_file, 2);

// some derivative tests
   dArg dummy(1.);
   std::vector<double> derivs;
   logLike.getFreeDerivs(dummy, derivs);
   
   std::vector<double> params;
   logLike.getFreeParamValues(params);
   
   double statval0 = logLike(dummy);
   
   for (unsigned int i = 0; i < params.size(); i++) {
      std::vector<double> new_params = params;
      double dparam = params[i]*1e-7;
      new_params[i] += dparam;
      logLike.setFreeParamValues(new_params);
      double statval = logLike(dummy);
      std::cout << derivs[i] << "  " 
                << (statval - statval0)/dparam << std::endl;
   }

// do the fit
   Lbfgs myOptimizer(logLike);
//   Minuit myOptimizer(logLike);

   int verbose = 3;
   try {
      myOptimizer.find_min(verbose);
   } catch (optimizers::Exception &eObj) {
      std::cerr << eObj.what() << std::endl;
   }

   std::vector<Parameter> parameters;
   logLike.getParams(parameters);

   for (unsigned int i = 0; i < parameters.size(); i++)
      std::cout << parameters[i].getName() << ": "
                << parameters[i].getValue() << std::endl;

   std::cout << "*** fit_anti_center: all tests completed ***\n" << std::endl;
}

void fit_3C279() {

   std::cout << "*** fit_3C279 ***" << std::endl;

   double ra = 193.98;
   double dec = -5.82;
   
// set the ROI cuts
   RoiCuts::setCuts(ra, dec, 30.);
   
/* read in the spacecraft data */
   std::string sc_file = test_path + "Data/one_src_sc_0000";
   int sc_hdu = 2;
   ScData::readData(sc_file, sc_hdu);
   
   logLike_ptsrc logLike;

   PointSource _3c279(ra, dec);
   
//    double Prefactor = 5.92e-5;
//    double Gamma = -1.96;
//    double Escale = 1;

   double Prefactor = 4.;
   double Gamma = -2.1;
   double Escale = 100.;
   
   PowerLaw pl(Prefactor, Gamma, Escale);

//set limits on index
   Parameter indexParam = pl.getParam("Index");
   indexParam.setBounds(-3.5, -1.);
   pl.setParam(indexParam);

// set limits on normalization
   Parameter prefactorParam = pl.getParam("Prefactor");
   prefactorParam.setBounds(1e-3, 1e3);  
   prefactorParam.setScale(1e-9);
   pl.setParam(prefactorParam);

   _3c279.setSpectrum(&pl);

   _3c279.setName("3C 279");

   logLike.addSource(&_3c279);

// read in the data
   std::string event_file = test_path + "Data/one_src_0000";
   logLike.getEvents(event_file, 2);

 // some derivative tests
   dArg dummy(1.);
   std::vector<double> derivs;
   logLike.getFreeDerivs(dummy, derivs);
   
   std::vector<double> params;
   logLike.getFreeParamValues(params);
   
   double statval0 = logLike(dummy);
   
   for (unsigned int i = 0; i < params.size(); i++) {
      std::vector<double> new_params = params;
      double dparam = params[i]*1e-7;
      new_params[i] += dparam;
      logLike.setFreeParamValues(new_params);
      double statval = logLike(dummy);
      std::cout << derivs[i] << "  " 
                << (statval - statval0)/dparam << std::endl;
   }

// do the fit using lbfgs_bcm
//   Lbfgs myOptimizer(logLike);
   Minuit myOptimizer(logLike);
   int verbose = 3;
   myOptimizer.find_min(verbose);

   std::vector<Parameter> parameters;
   logLike.getParams(parameters);
   
   for (unsigned int i = 0; i < parameters.size(); i++)
      std::cout << parameters[i].getName() << ": "
                << parameters[i].getValue() << std::endl;

   std::cout << "*** fit_3C279: all tests completed ***\n" << std::endl;
}

/***********************/
/* logLike_ptsrc tests */
/***********************/
void test_logLike_ptsrc() {

   std::cout << "*** test_logLike_ptsrc ***" << std::endl;

   double ra = 193.98;
   double dec = -5.82;

// set the ROI cuts
   RoiCuts::setCuts(ra, dec, 30.);

/* read in the spacecraft data */
   std::string sc_file = test_path + "Data/one_src_sc_0000";
   int sc_hdu = 2;
   ScData::readData(sc_file, sc_hdu);

   logLike_ptsrc logLike;

   PointSource _3c279(ra, dec);

   double Prefactor = 7.;
   double Gamma = -1.96;
   double Escale = 100.;

   PowerLaw pl(Prefactor, Gamma, Escale);
   Parameter prefactorParam = pl.getParam("Prefactor");
   prefactorParam.setScale(1e-9);
   pl.setParam(prefactorParam);

   _3c279.setSpectrum(&pl);

   _3c279.setName("3C 279");

   logLike.addSource(&_3c279);

   std::string event_file = test_path + "Data/one_src_0000";
   logLike.getEvents(event_file, 2);

   report_SrcModel_values(logLike);

   std::vector<Parameter> params;
   logLike.getFreeParams(params);

// vary over the Prefactor of the Spectrum and compute the
// log-likelihood at each step

   double xmin = 3.;
   double xmax = 10.;
   int nx = 20;
   double xstep = log(xmax/xmin)/(nx-1.);

   unsigned int j;
   for (int i = 0; i < nx; i++) {
      for (j = 0; j < params.size(); j++) {
         if (params[j].getName() == "Prefactor") {
            params[j].setValue(xmin*exp(xstep*i));
            break;
         }
      }
      std::vector<double> param_vals;
      for (unsigned int k = 0; k < params.size(); k++)
         param_vals.push_back(params[k].getValue());
      logLike.setFreeParamValues(param_vals);
      dArg dummy(1.);
      std::cout << params[j].getValue() << "  " 
                << logLike(dummy) << std::endl;
   }

   std::cout << "*** test_logLike_ptsrc: all tests completed ***\n" 
             << std::endl;

}
// logLike_ptsrc tests

/********************/
/* Psf class tests */
/********************/
void test_Psf_class() {

   std::cout << "*** test_Psf_class ***" << std::endl;

   Psf *psf = Psf::instance();

   int nenergy = 10;
   double emin = 30;
   double emax = 3e4;
   double estep = log(emax/emin)/(nenergy - 1);

   std::cout << "\nPoint Spread Function data:" << std::endl;
   for (double inc = 0.; inc <= 70.; inc += 10.) {
      for (int i = 0; i < nenergy; i++) {
         double energy = emin*exp(estep*i);
         std::cout << energy << "  ";

         std::vector<double> psf_params;
         psf->fillPsfParams(energy, inc, psf_params);
         std::cout << psf_params[0] << "  "
                   << psf_params[1] << "  "
                   << psf_params[2] << "\n";
      }
      std::cout << std::endl;
   }
   std::cout << "*** test_Psf_class: all tests completed ***\n" << std::endl;

} // Psf class tests

/********************/
/* Aeff class tests */
/********************/
void test_Aeff_class() {

   std::cout << "*** test_Aeff_class ***" << std::endl;

   Aeff *aeff = Aeff::instance();

   int nenergy = 10;
   double emin = 30;
   double emax = 3e4;
   double estep = log(emax/emin)/(nenergy - 1);

   std::cout << "\nEffective area data:" << std::endl;
   for (int i = 0; i < nenergy; i++) {
      double energy = emin*exp(estep*i);
      std::cout << energy << "  ";
      for (double inc = 0.; inc <= 70.; inc += 10.) {
         std::cout << (*aeff)(energy, inc) << "  ";
      }
      std::cout << std::endl;
   }
   std::cout << std::endl;

   std::cout << "*** test_Aeff_class: all tests completed ***\n" << std::endl;

} // Aeff class tests

/***************************/
/* PointSource class tests */
/***************************/
void test_PointSource_class() {

   std::cout << "*** test_PointSource_class ***" << std::endl;

// set the ROI cuts
   RoiCuts::setCuts(193.98, -5.82, 30.);

/* read in the spacecraft data */
   std::string sc_file = test_path + "Data/one_src_sc_0000";
   int sc_hdu = 2;
   ScData::readData(sc_file, sc_hdu);

/* put this source at the center of the extraction region for
   Data/one_src_0000 (ra, dec) = (193.98, -5.82) */
   
   PointSource my_ptsrc(193.98, -5.82);

   my_ptsrc.setName("3C_279");

/* create a PowerLaw object for the source spectrum */

   PowerLaw source_pl(1., -2.1, 100.); // f(E) = (E/0.1 GeV)^(-2.1)
   my_ptsrc.setSpectrum(&source_pl);

/* compute the source flux at various energies and aspect angles */

   int nenergy = 10;
   double emin = 30.;
   double emax = 3e4;
   double estep = log(emax/emin)/(nenergy - 1);

   std::cout << "\nThe on-source counts spectrum of "
             << my_ptsrc.getName() 
             << " as a function of time: " 
             << std::endl;

   double photon_flux = 0;

   double maxtime = 95.*60.;
   double tstep = 1e2;
   for (double time = 0; time < maxtime; time += tstep) {
      std::cout << "Time = " << time << std::endl;
      for (int i = 0; i < nenergy; i++) { // loop over energies
         double energy = emin*exp(estep*i);
         if (i == 0) photon_flux = 
                        my_ptsrc.fluxDensity(energy, time, my_ptsrc.getDir());
         if (photon_flux > 0) {
            std::cout << energy << "  ";
            std::cout << my_ptsrc.fluxDensity(energy, time, my_ptsrc.getDir())
                      << std::endl;
         }
      }
      if (photon_flux > 0) std::cout << std::endl;
   }

   std::cout << "*** test_PointSource_class: all tests completed ***\n" 
             << std::endl;

} // PointSource class tests

/*********************/
/* Event class tests */
/*********************/
void test_Event_class() {

   std::cout << "*** test_Event_class ***" << std::endl;

// read in the event data, then stuff them into an Event class vector

   Table evt_table;

/* read in EVENT file */
   std::string event_file = test_path + "Data/one_src_0000";

   evt_table.add_columns("RA DEC energy time SC_x SC_y SC_z zenith_angle");
   evt_table.read_FITS_table(event_file, 2);

   typedef std::pair<long, double*> tableColumn;
   tableColumn ra(evt_table[0].dim, evt_table[0].val);
   tableColumn dec(evt_table[1].dim, evt_table[1].val);
   tableColumn energy(evt_table[2].dim, evt_table[2].val);
   tableColumn time(evt_table[3].dim, evt_table[3].val);
   tableColumn sc_x(evt_table[4].dim, evt_table[4].val);
   tableColumn sc_y(evt_table[5].dim, evt_table[5].val);
   tableColumn sc_z(evt_table[6].dim, evt_table[6].val);
   tableColumn zenangle(evt_table[7].dim, evt_table[7].val);

   std::vector<Event> my_events;

// get pointer to RoiCuts
   RoiCuts *roi_cuts = RoiCuts::instance();
   RoiCuts::setCuts(193.98, -5.82, 30);

   unsigned int nReject = 0;

   for (int i = 0; i < ra.first; i++) {
// compute sc_ra and sc_dec from direction cosines (it would be nice
// if SkyDir could do this for me....)
      double sc_ra = atan2(sc_y.second[i], sc_x.second[i])*180./M_PI;
      double sc_dec = asin(sc_z.second[i])*180./M_PI;

      Event thisEvent(ra.second[i], dec.second[i], energy.second[i],
                      time.second[i], sc_ra, sc_dec,
                      cos(zenangle.second[i]*M_PI/180.));
      if (roi_cuts->accept(thisEvent)) {
         my_events.push_back(thisEvent);
      } else {
         nReject++;
      }
   }

   std::cout << "Out of " << ra.first << " events, "
             << my_events.size() << " were accepted, and "
             << nReject << " were rejected.\n" << std::endl;

   astro::SkyDir _3c279_dir(193.98, -5.82);

/* write out what's available for nmax of these guys */
// since I know these photons were simulated for 3C 279, write out the
// difference for later binning to get the psf

   int nmax = 20;
   for (int i = 0; i < nmax; i++) {
      astro::SkyDir thisEventDir = my_events[i].getDir();
      std::cout << thisEventDir.ra() << "  " 
                << thisEventDir.dec() << "  "
                << my_events[i].getEnergy() << "  "
                << my_events[i].getArrTime() << "  "
                << my_events[i].getSeparation(_3c279_dir)*180./M_PI
                << std::endl;
   }
   std::cout << std::endl;

   std::cout << "*** test_Event_class: all tests completed ***\n" << std::endl;

} // Event class tests

/*********************/
/* Table class tests */
/*********************/
void test_Table_class() {

   std::cout << "*** test_Table_class ***" << std::endl;

/* read in PSF parameters */
   std::string psf_file = test_path + "CALDB/psf_lat.fits";

   Table psf_data;
      
   psf_data.add_columns("ENERGY THETA SIG1_F SIG2_F W");
   psf_data.read_FITS_table(psf_file, 2);

   int nenergy = psf_data[0].dim;
   double *energy = psf_data[0].val;
   int ntheta = psf_data[1].dim;
   double *theta = psf_data[1].val;

   std::cout << "\nFrom " << psf_file << ": \n";

   std::cout << "energies: ";
   for (int i = 0; i < nenergy; i++) 
      std::cout << energy[i] << "  ";
   std::cout << std::endl;

   std::cout << "theta values: ";
   for (int i = 0; i < ntheta; i++) 
      std::cout << theta[i] << "  ";
   std::cout << std::endl;
   std::cout << std::endl;

/* read in AEFF parameters */
   std::string aeff_file = test_path + "CALDB/aeff_lat.fits";

   Table aeff_data;

   aeff_data.add_columns("ENERGY THETA AEFF_F");
   aeff_data.read_FITS_table(aeff_file, 2);

   nenergy = aeff_data[0].dim;
   energy = aeff_data[0].val;
   ntheta = aeff_data[1].dim;
   theta = aeff_data[1].val;
   
   std::cout << "From " << aeff_file << ": \n";
   std::cout << "energies: ";
   for (int i = 0; i < nenergy; i++) 
      std::cout << energy[i] << "  ";
   std::cout << std::endl;

   std::cout << "theta values: ";
   for (int i = 0; i < ntheta; i++) 
      std::cout << theta[i] << "  ";
   std::cout << std::endl;
   std::cout << std::endl;

/* read in EVENT file */
   std::string event_file = test_path + "Data/one_src_0000";

   Table event_data;

   event_data.add_columns("RA DEC energy time SC_x SC_y SC_z zenith_angle");
   event_data.read_FITS_table(event_file, 2);

   int nRA = event_data[0].dim;
   double* RA = event_data[0].val;
   int nDEC = event_data[1].dim;
   double* DEC = event_data[1].val;
   nenergy = event_data[2].dim;
   energy = event_data[2].val;
   int ntime = event_data[3].dim;
   double* time = event_data[3].val;

   int nmax = 5;

   std::cout << "Some of the data from " << event_file << ":\n";
   std::cout << nRA << "  " << nDEC << "  "
             << nenergy << "  " << ntime << std::endl;
   for (int i = 0; i < nmax; i++) {
      std::cout << RA[i] << "  "
                << DEC[i] << "  "
                << energy[i] << "  "
                << time[i] << std::endl;
   }
   std::cout << std::endl;

/* Spacecraft data file */
   std::string sc_file = test_path + "Data/one_src_sc_0000";

   Table sc_data;

   sc_data.add_columns("time SC_x SC_y SC_z SAA_flag");
   sc_data.read_FITS_table(sc_file, 2);

   ntime= sc_data[0].dim;
   time = sc_data[0].val;
   int nSC_x = sc_data[1].dim;
   double* SC_x = sc_data[1].val;
   int nSC_y = sc_data[2].dim;
   double* SC_y = sc_data[2].val;
   int nSC_z = sc_data[3].dim;
   double* SC_z = sc_data[3].val;

   std::cout << "Some of the data from " << sc_file << ":\n";
   std::cout << ntime << "  " << nSC_x << "  "
             << nSC_y << "  " << nSC_z << std::endl;
   for (int i = 0; i < nmax; i++) {
      std::cout << time[i] << "  "
                << SC_x[i] << "  "
                << SC_y[i] << "  "
                << SC_z[i] << std::endl;
   }
   std::cout << std::endl;

   std::cout << "*** test_Table_class: all tests completed ***\n" << std::endl;

} // Table class tests

/***************************/
/* SourceModel class tests */
/***************************/
void test_SourceModel_class() {

   std::cout << "*** test_SourceModel_class ***" << std::endl;
   
   SourceModel SrcModel;
   bool computeExposure = false;
   
/* create some point sources */
   PointSource _3c279;
   _3c279.setDir(193.98, -5.82, computeExposure);
   _3c279.setSpectrum(new PowerLaw(74.2, -1.96, 0.1));
   _3c279.setName("3C 279");

   PointSource _3c273;
   _3c273.setDir(187.25, 2.17, computeExposure);
   _3c273.setSpectrum(new PowerLaw(15.4, -2.58, 0.1));
   _3c273.setName("3C 273");

   PointSource Crab;
   Crab.setDir(83.57, 22.01, computeExposure);
   Crab.setSpectrum(new PowerLaw(226.2, -2.19, 0.1));
   Crab.setName("Crab Pulsar");

   PointSource Vela;
   Vela.setDir(128.73, -45.20, computeExposure);
   Vela.setSpectrum(new PowerLaw(834.3, -1.69, 0.1));
   Vela.setName("Vela Pulsar");

/* add these guys to the SouceModel */
   SrcModel.addSource(&_3c279);
   SrcModel.addSource(&_3c273);
   SrcModel.addSource(&Crab);
   SrcModel.addSource(&Vela);
   report_SrcModel_values(SrcModel);

/* delete a Source */
   SrcModel.deleteSource("Crab Pulsar");
   report_SrcModel_values(SrcModel);

/* add another function */
   PointSource Geminga;
   Geminga.setDir(98.49, 17.86, computeExposure);
   Geminga.setSpectrum(new PowerLaw(352.9, -1.66, 0.1));
   Geminga.setName("Geminga");
   SrcModel.addSource(&Geminga);
   
/* make its scale factor free (this could be made easier, e.g., by
   giving direct parameter access from PointSource) */
   Parameter param = SrcModel.getParam("Scale", "Spectrum", "Geminga");
   param.setFree(true);
   SrcModel.setParam(param, "Spectrum", "Geminga");
   report_SrcModel_values(SrcModel);

/* derivative tests */
   dArg x(20.);
   std::vector<double> params_save;
   SrcModel.getParamValues(params_save);
   std::vector<double> sm_derivs;
   SrcModel.getDerivs(x, sm_derivs);
   std::vector<double> params = params_save;

   std::cout << "SrcModel Derivatives: " << std::endl;
   for (unsigned int i = 0; i < sm_derivs.size(); i++) {
      std::cout << sm_derivs[i] << ":  ";

// compute the numerical derivative wrt this parameter
      double delta_param = fabs(params_save[i]/1e7);
      double num_deriv = SrcModel(x);
      params[i] += delta_param;
      SrcModel.setParamValues(params);
      num_deriv -= SrcModel(x);
      num_deriv /= delta_param;
      std::cout << -num_deriv << std::endl;

// reset the parameters for next time around
      SrcModel.setParamValues(params_save);
      params = params_save;
   }
   std::cout << std::endl;

/* derivative tests for free parameters only */
   SrcModel.getFreeParamValues(params_save);
   SrcModel.getFreeDerivs(x, sm_derivs);
   params = params_save;

   std::cout << "SrcModel Derivatives for Free Parameters: " << std::endl;
   for (unsigned int i = 0; i < sm_derivs.size(); i++) {
      std::cout << sm_derivs[i] << ":  ";

// compute the numerical derivative wrt this parameter
      double delta_param = fabs(params_save[i]/1e7);
      double num_deriv = SrcModel(x);
      params[i] += delta_param;
      SrcModel.setFreeParamValues(params);
      num_deriv -= SrcModel(x);
      num_deriv /= delta_param;
      std::cout << -num_deriv << std::endl;

// reset the parameters for next time around
      SrcModel.setFreeParamValues(params_save);
      params = params_save;
   }
   std::cout << std::endl;

   std::cout << "*** test_SourceModel_class: all tests completed ***\n" 
             << std::endl;

} // Source Model Class tests

void report_SrcModel_values(const SourceModel &SrcModel) {

/* retrieve Free parameter values */
  std::vector<double> paramValues;
  SrcModel.getFreeParamValues(paramValues);
  for (unsigned int i = 0; i < paramValues.size(); i++) {
     std::cout << paramValues[i] << "  ";
  }
  std::cout << std::endl;

/* retrieve all parameter values */
  SrcModel.getParamValues(paramValues);
  for (unsigned int i = 0; i < paramValues.size(); i++) {
     std::cout << paramValues[i] << "  ";
  }
  std::cout << std::endl;
   
/* retrieve Source names */
  std::vector<std::string> srcNames;
  SrcModel.getSrcNames(srcNames);
  for (unsigned int i = 0; i < srcNames.size(); i++) {
     std::cout << srcNames[i] << "  ";
  }
  std::cout << std::endl;
  std::cout << std::endl;
}

void read_SC_Response_data() {

/* get root path to test data */   
   const char * root = ::getenv("LIKELIHOODROOT");
   if (!root) {  //use relative path from cmt directory
      root_path = "..";
   } else {
      root_path = std::string(root);
   }
   test_path = root_path + "/src/test/";

/* instantiate the Psf and read in its data */
   Psf * psf = Psf::instance();
   std::string psf_file = test_path + "CALDB/psf_lat.fits";
   try {
      psf->readPsfData(psf_file, Response::Combined);
   } catch (Exception &eObj) {
      std::cerr << eObj.what() << std::endl;
      throw;
   }

/* instantiate the Aeff and read in its data */
   Aeff * aeff = Aeff::instance();
   std::string aeff_file = test_path + "CALDB/aeff_lat.fits";
   aeff->readAeffData(aeff_file, Response::Combined);
}

void print_fit_results(SourceModel &stat) {
   std::vector<std::string> srcNames;
   stat.getSrcNames(srcNames);
   std::vector<Parameter> parameters;
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
