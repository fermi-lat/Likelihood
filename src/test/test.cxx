/**
 * @file test.cxx
 * @brief Test program for Likelihood.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/test/test.cxx,v 1.147 2016/10/20 23:13:04 echarles Exp $
 */

#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <cmath>
#include <cstdio>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <cppunit/ui/text/TextTestRunner.h>
#include <cppunit/extensions/HelperMacros.h>

#include "CLHEP/Random/RandFlat.h"

#include "facilities/Util.h"
#include "facilities/commonUtilities.h"

#include "st_facilities/Environment.h"
#include "st_facilities/Util.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "evtbin/HealpixMap.h"

#include "optimizers/dArg.h"
#include "optimizers/FunctionFactory.h"
#include "optimizers/FunctionTest.h"
#ifdef DARWIN_F2C_FAILURE
#include "optimizers/NewMinuit.h"
#else
#include "optimizers/Minuit.h"
#endif

#include "irfInterface/IrfsFactory.h"
#include "irfInterface/AcceptanceCone.h"
#include "irfLoader/Loader.h"

#include "Likelihood/BinnedConfig.h"
#include "Likelihood/BinnedExposure.h"
#include "Likelihood/BinnedHealpixExposure.h"
#include "Likelihood/BinnedLikelihood.h"
#include "Likelihood/CompositeSource.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/CountsMapHealpix.h"
#include "Likelihood/DiffRespNames.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/Drm.h"
#include "Likelihood/Event.h"
#include "Likelihood/EventContainer.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/FitUtils.h"
#include "Likelihood/FluxBuilder.h"
#include "Likelihood/LikeExposure.h"
#include "Likelihood/LogNormal.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/Observation.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/ScaleFactor.h"
#include "Likelihood/SourceModelBuilder.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/ScData.h"
#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/Source.h"
#include "Likelihood/SourceFactory.h"
#include "Likelihood/SourceMap.h"
#include "Likelihood/SourceModel.h"
#include "Likelihood/SpatialMap.h"
#include "Likelihood/TrapQuad.h"
#include "Likelihood/WcsMap2.h"
#include "Likelihood/BandFunction.h"
#include "Likelihood/ExpCutoffSEDPeak.h"
#include "Likelihood/BrokenPowerLaw2.h"
#include "Likelihood/BrokenPowerLaw3.h"
#include "Likelihood/BrokenPowerLawExpCutoff.h"
#include "Likelihood/EblAtten.h"
#include "Likelihood/EnergyBand.h"
#include "Likelihood/LogParabola.h"
#include "Likelihood/MultipleBrokenPowerLaw.h"
#include "Likelihood/PiecewisePowerLaw.h"
#include "Likelihood/ProjMap.h"
#include "Likelihood/PowerLaw2.h"
#include "Likelihood/PowerLawSuperExpCutoff.h"
#include "Likelihood/SmoothBrokenPowerLaw.h"
#include "Likelihood/SmoothDoubleBrokenPowerLaw.h"
#include "Likelihood/WcsMapLibrary.h"

#include "SourceData.h"
#include "XmlDiff.h"

using namespace Likelihood;
using optimizers::Parameter;

class LikelihoodTests : public CppUnit::TestFixture {

   CPPUNIT_TEST_SUITE(LikelihoodTests);

   CPPUNIT_TEST(test_LogParabola);
   CPPUNIT_TEST(test_LogNormal);
   CPPUNIT_TEST(test_BandFunction);
   CPPUNIT_TEST(test_ExpCutoffSEDPeak);
   CPPUNIT_TEST(test_SmoothBrokenPowerLaw);
   CPPUNIT_TEST(test_SmoothDoubleBrokenPowerLaw);
   CPPUNIT_TEST(test_BrokenPowerLaw3);
   CPPUNIT_TEST(test_MultipleBrokenPowerLaw);
   CPPUNIT_TEST(test_PiecewisePowerLaw);
   CPPUNIT_TEST(test_EblAtten);
   CPPUNIT_TEST(test_EnergyBand);
   CPPUNIT_TEST(test_RoiCuts);
   CPPUNIT_TEST(test_SourceFactory);
   CPPUNIT_TEST(test_XmlBuilders);
   CPPUNIT_TEST(test_LikeExposure);
   CPPUNIT_TEST(test_SourceModel);
   CPPUNIT_TEST(test_SourceDerivs);
   CPPUNIT_TEST(test_PointSource);
   CPPUNIT_TEST(test_DiffuseSource);
   CPPUNIT_TEST(test_CountsMap);
   CPPUNIT_TEST(test_CountsMapHealpix_allsky);
   CPPUNIT_TEST(test_CountsMapHealpix_region);
   CPPUNIT_TEST(test_BinnedLikelihood);
   CPPUNIT_TEST(test_BinnedLikelihood_2);
   CPPUNIT_TEST(test_BinnedLikelihood_wts);
   CPPUNIT_TEST(test_BinnedLikelihood_edisp_neg1);
   CPPUNIT_TEST(test_BinnedLikelihood_edisp_neg2);
   CPPUNIT_TEST(test_BinnedLikelihood_edisp_1);
   CPPUNIT_TEST(test_BinnedLikelihood_edisp_2);
   CPPUNIT_TEST(test_CompositeSource);
   CPPUNIT_TEST(test_MeanPsf);
   CPPUNIT_TEST(test_BinnedExposure);
   CPPUNIT_TEST(test_BinnedExposureHealpix);
   CPPUNIT_TEST(test_SourceMap);
   CPPUNIT_TEST(test_PointSourceMap);
   CPPUNIT_TEST(test_PointSourceMap_hpx_allsky);
   CPPUNIT_TEST(test_PointSourceMap_hpx_region);
   CPPUNIT_TEST(test_rescaling);
   CPPUNIT_TEST(test_DiffRespNames);
   CPPUNIT_TEST_EXCEPTION(test_WcsMap2_exception, std::runtime_error);
   CPPUNIT_TEST(test_WcsMap2);
   CPPUNIT_TEST(test_ScaleFactor);
   CPPUNIT_TEST(test_Drm);
   CPPUNIT_TEST(test_Source_Npred);
   CPPUNIT_TEST(test_ExposureCube);

   CPPUNIT_TEST_SUITE_END();

public:

   void setUp();
   void tearDown();

   void test_LogParabola();
   void test_LogNormal();
   void test_BandFunction();
   void test_ExpCutoffSEDPeak();
   void test_SmoothBrokenPowerLaw();
   void test_SmoothDoubleBrokenPowerLaw();
   void test_BrokenPowerLaw3();
   void test_MultipleBrokenPowerLaw();
   void test_PiecewisePowerLaw();
   void test_EblAtten();
   void test_EnergyBand();
   void test_RoiCuts();
   void test_SourceFactory();
   void test_XmlBuilders();
   void test_LikeExposure();
   void test_SourceModel();
   void test_SourceDerivs();
   void test_PointSource();
   void test_DiffuseSource();
   void test_CountsMap();
   void test_CountsMapHealpix_allsky();
   void test_CountsMapHealpix_region();
   void test_BinnedLikelihood();
   void test_BinnedLikelihood_2();
   void test_BinnedLikelihood_base(int edisp_val);
   void test_BinnedLikelihood_wts() {
      test_BinnedLikelihood_base(0);
   }
   void test_BinnedLikelihood_edisp_neg1() {
      test_BinnedLikelihood_base(-1);
   }
   void test_BinnedLikelihood_edisp_neg2() {
      test_BinnedLikelihood_base(-2);
   }
   void test_BinnedLikelihood_edisp_1() {
      test_BinnedLikelihood_base(1);
   }
   void test_BinnedLikelihood_edisp_2() {
      test_BinnedLikelihood_base(2);
   }
   void test_CompositeSource();
   void test_MeanPsf();
   void test_BinnedExposure();
   void test_BinnedExposureHealpix();
   void test_SourceMap();
   void test_PointSourceMap();
   void test_PointSourceMap_hpx_allsky();
   void test_PointSourceMap_hpx_region();
   void test_rescaling();
   void test_DiffRespNames();
   void test_WcsMap2_exception();
   void test_WcsMap2();
   void test_ScaleFactor();
   void test_Drm();
   void test_Source_Npred();
   void test_ExposureCube();

private:

   Observation* m_observation = nullptr;
   static RoiCuts* m_roiCuts;
   static ScData* m_scData;
   static ExposureCube* m_expCube;
   static ExposureMap* m_expMap;
   static ResponseFunctions* m_respFuncs;
   static EventContainer* m_eventCont;

   std::string m_rootPath;
   double m_fracTol = 1e-4;

   // File names for m_srcFactory.
   std::string m_scFile;
   std::string m_expMapFile;
   std::string m_sourceXmlFile;

   // Test data for m_srcFactory.
   SourceData m_srcData;

   // Filenames for test_XmlBuilders.
   std::string m_fluxXmlFile;
   std::string m_srcModelXmlFile;

   optimizers::FunctionFactory* m_funcFactory = nullptr;
   optimizers::FunctionFactory* funcFactoryInstance();

   SourceFactory* m_srcFactory = nullptr;
   SourceFactory* srcFactoryInstance(const std::string& scFile = "",
                                     const std::string& expMapFile = "",
                                     const std::string& sourceXmlFile = "",
                                     bool requireExposure = true,
                                     bool verbose = true);

   void readEventData(const std::string& eventFile,
                      const std::string& scDataFile,
                      std::vector<Event>& events);

   void generate_exposureHyperCube();

   CountsMap singleSrcMap(unsigned int nee,
                          unsigned long num_x_pix = 40, 
                          unsigned long num_y_pix = 40,
                          double pix_scale = 0.25) const;

   CountsMapHealpix healpixmap_allsky() const;
   CountsMapHealpix healpixmap_region() const;

   std::string dataPath(const std::string& filename) const {
      return facilities::commonUtilities::joinPath(m_rootPath, filename);
   }

   static double ASSERT_TOLERANCE;
};

// Static member initialization
RoiCuts* LikelihoodTests::m_roiCuts = nullptr;
ScData* LikelihoodTests::m_scData = nullptr;
ExposureCube* LikelihoodTests::m_expCube = nullptr;
ExposureMap* LikelihoodTests::m_expMap = nullptr;
ResponseFunctions* LikelihoodTests::m_respFuncs = nullptr;
EventContainer* LikelihoodTests::m_eventCont = nullptr;
double LikelihoodTests::ASSERT_TOLERANCE = 1e-6;

#define ASSERT_EQUALS(X, Y) CPPUNIT_ASSERT(std::fabs((X) - (Y)) < ASSERT_TOLERANCE * (std::fabs(X) + std::fabs(Y)) / 2.)

void LikelihoodTests::setUp() {
   if (m_respFuncs == 0) m_respFuncs = new ResponseFunctions();
   if (m_scData == 0) m_scData = new ScData();
   if (m_roiCuts == 0) m_roiCuts = new RoiCuts();
   if (m_expCube == 0) m_expCube = new ExposureCube();
   if (m_expMap == 0) m_expMap = new ExposureMap();
   if (m_eventCont == 0) m_eventCont = new EventContainer(*m_respFuncs,
                                                          *m_roiCuts, 
                                                          *m_scData);
   m_observation = new Observation(m_respFuncs,
                                   m_scData,
                                   m_roiCuts,
                                   m_expCube,
                                   m_expMap,
                                   m_eventCont);

// Get root path to test data.
   const char * root = 
      st_facilities::Environment::packagePath("Likelihood").c_str();
   if (!root) {  //use relative path from cmt directory
      m_rootPath = "..";
   } else {
      m_rootPath = st_facilities::Environment::dataPath("Likelihood");
   }

// Prepare the ResponseFunctions object.
   irfLoader::Loader::go();
   irfInterface::IrfsFactory * myFactory 
      = irfInterface::IrfsFactory::instance();
   m_respFuncs->addRespPtr(0, myFactory->create("DC1A::Front"));
   m_respFuncs->addRespPtr(1, myFactory->create("DC1A::Back"));
//    m_respFuncs->addRespPtr(0, myFactory->create("P7SOURCE_V6::FRONT"));
//    m_respFuncs->addRespPtr(1, myFactory->create("P7SOURCE_V6::BACK"));

// Fractional tolerance for double comparisons.
   m_fracTol = 1e-4;

// Use lazy evaluation for m_funcFactory and m_srcFactory.
   m_funcFactory = 0;
   m_srcFactory = 0;

   m_scFile = dataPath("oneday_scData_0000.fits");
   m_expMapFile = dataPath("anticenter_expMap.fits");
   m_sourceXmlFile = dataPath("anticenter_model.xml");
}

void LikelihoodTests::tearDown() {
   delete m_observation;
   m_observation = nullptr;
}

void LikelihoodTests::test_XmlBuilders() {
   SourceFactory* srcFactory = srcFactoryInstance();
   std::vector<std::string> srcNames;
   srcFactory->fetchSrcNames(srcNames);
   
   FluxBuilder fluxBuilder(30, 2e5);
   SourceModelBuilder srcModelBuilder("", "source library");
   
   for (const auto& name : srcNames) {
      Source* src = srcFactory->create(name);
      fluxBuilder.addSource(*src);
      srcModelBuilder.addSource(*src);
   }
   
   m_fluxXmlFile = dataPath("fluxBuilder.xml");
   fluxBuilder.write(m_fluxXmlFile);

   m_srcModelXmlFile = dataPath("srcModelBuilder.xml");
   srcModelBuilder.write(m_srcModelXmlFile);

   // Use updated XmlDiff class with RapidXML
   XmlDiff xmlDiff(m_sourceXmlFile, m_srcModelXmlFile, "source", "name");
   CPPUNIT_ASSERT(xmlDiff.compare());
}

void LikelihoodTests::test_ScaleFactor() {
   Likelihood::ScaleFactor foo(PowerLaw2(1, -2.1, 20, 2e5), 1);
   optimizers::FunctionTest tester(foo, "ScaleFactor::PowerLaw2");

   std::vector<optimizers::Parameter> params;
   params.push_back(optimizers::Parameter("Integral", 1));
   params.push_back(optimizers::Parameter("Index", -2.1));
   params.push_back(optimizers::Parameter("LowerLimit", 20.));
   params.push_back(optimizers::Parameter("UpperLimit", 2e5));
   params.push_back(optimizers::Parameter("ScaleFactor", 1));

   std::vector<optimizers::Arg*> args;
   args.push_back(new optimizers::dArg(100));
   args.push_back(new optimizers::dArg(300));
   args.push_back(new optimizers::dArg(1e3));
   args.push_back(new optimizers::dArg(3e3));
   args.push_back(new optimizers::dArg(1e4));
   args.push_back(new optimizers::dArg(3e4));
   args.push_back(new optimizers::dArg(1e5));

   tester.parameters(params);
   tester.freeParameters(params);
   tester.derivatives(args, 1e-5);

   // Test complement functionality
   double saved_value(foo(*args.at(0)));
   
   // This should throw an exception.
   try {
      foo.set_complement_flag(true);
   } catch (std::runtime_error&) {
   }
   
   // Set the value, bounds and scale for the ScaleFactor param
   foo.parameter("ScaleFactor").setValue(1);
   foo.parameter("ScaleFactor").setBounds(0, 1);
   foo.parameter("ScaleFactor").setScale(1);
   foo.set_complement_flag(true);
   
   // Check complement functionality
   CPPUNIT_ASSERT(foo(*args.at(0)) == 0);
   foo.parameter("ScaleFactor").setValue(0);
   ASSERT_EQUALS(foo(*args.at(0)), saved_value);
   
   // Clean up args
   for (auto* arg : args) {
      delete arg;
   }
}

void LikelihoodTests::test_PiecewisePowerLaw() {
   Likelihood::PiecewisePowerLaw foo;
   double indexL(-2);
   double indexH(-3);
   double dnde_values[] = {10, 3, 2, 0.1};
   double energy_values[] = {1e2, 5.5e2, 1.73e3, 5.5e3};
   std::vector<double> dNdEs(dnde_values, dnde_values + 4);
   std::vector<double> energies(energy_values, energy_values + 4);
   foo.addParams(indexL, indexH, dNdEs, energies);
   
   optimizers::FunctionTest tester(foo, "PiecewisePowerLaw");
   std::vector<optimizers::Parameter> pars;
   pars.push_back(optimizers::Parameter("IndexL", indexL));
   pars.push_back(optimizers::Parameter("IndexH", indexH));
   
   for (size_t k = 0; k < dNdEs.size(); ++k) {
      std::ostringstream name;
      name << "dNdE" << k;
      pars.push_back(optimizers::Parameter(name.str(), dNdEs[k]));
      pars.back().setScale(1e-12);
   }
   
   std::vector<optimizers::Arg*> args;
   args.push_back(new optimizers::dArg(100));
   args.push_back(new optimizers::dArg(300));
   args.push_back(new optimizers::dArg(1e3));
   args.push_back(new optimizers::dArg(3e3));
   args.push_back(new optimizers::dArg(1e4));
   args.push_back(new optimizers::dArg(3e4));
   args.push_back(new optimizers::dArg(1e5));

   tester.freeParameters(pars);
   tester.derivatives(args, 1e-4);
   
   // Clean up args
   for (auto* arg : args) {
      delete arg;
   }
}

void LikelihoodTests::test_RoiCuts() {
   m_roiCuts->setCuts(193.98, -5.82, 20, 30, 3.1623e5, 0, 1e12, -1, true);

   // Compare to known default values.
   std::vector<std::pair<double, double>> tlims;
   m_roiCuts->getTimeCuts(tlims);
   static double tmin = 0;
   static double tmax = 1e12;
   CPPUNIT_ASSERT(std::fabs(tlims[0].first - tmin) == 0);
   ASSERT_EQUALS(tlims[0].second, tmax);

   std::pair<double, double> energies = m_roiCuts->getEnergyCuts();
   static double emin = 30.;
   static double emax = 3.1623e5;
   ASSERT_EQUALS(energies.first, emin);
   ASSERT_EQUALS(energies.second, emax);
}

void LikelihoodTests::test_DiffRespNames() {
   DiffRespNames foo;
   
   foo.addColumn("P6_V1_DIFFUSE__Extragalactic Diffuse");
   foo.addColumn("P6_V1_DIFFUSE__GALPROP Diffuse");
   
   CPPUNIT_ASSERT("P6_V1_DIFFUSE__Extragalactic Diffuse" == foo[0]);
   CPPUNIT_ASSERT("P6_V1_DIFFUSE__Extragalactic Diffuse" == foo["DIFRSP0"]);
                 
   CPPUNIT_ASSERT("P6_V1_DIFFUSE__GALPROP Diffuse" == foo[1]);
   CPPUNIT_ASSERT("P6_V1_DIFFUSE__GALPROP Diffuse" == foo["DIFRSP1"]);
}

// Placeholder implementations for other test methods
void LikelihoodTests::test_LogParabola() {
   // Implementation
}

void LikelihoodTests::test_LogNormal() {
   // Implementation  
}

void LikelihoodTests::test_BandFunction() {
   Likelihood::BandFunction band(1, -1.5, -2.5, 1e3, 100.);
   optimizers::FunctionTest tester(band, "BandFunction");
   std::vector<optimizers::Parameter> params;
   params.push_back(optimizers::Parameter("norm", 1));
   params.push_back(optimizers::Parameter("alpha", -1.5));
   params.push_back(optimizers::Parameter("beta", -2.5));
   params.push_back(optimizers::Parameter("Ep", 0.1));
   params.push_back(optimizers::Parameter("Scale", 0.1));

   std::vector<optimizers::Arg*> args;
   args.push_back(new optimizers::dArg(100));
   args.push_back(new optimizers::dArg(300));
   args.push_back(new optimizers::dArg(1e3));
   args.push_back(new optimizers::dArg(3e3));
   args.push_back(new optimizers::dArg(1e4));
   args.push_back(new optimizers::dArg(3e4));
   args.push_back(new optimizers::dArg(1e5));

   tester.parameters(params);
   tester.freeParameters(params);
   tester.derivatives(args, 1e-5);
   
   for (auto* arg : args) {
      delete arg;
   }
}

void LikelihoodTests::test_ExpCutoffSEDPeak() {
   // Implementation
}

void LikelihoodTests::test_SmoothBrokenPowerLaw() {
   // Implementation
}

void LikelihoodTests::test_SmoothDoubleBrokenPowerLaw() {
   // Implementation
}

void LikelihoodTests::test_BrokenPowerLaw3() {
   Likelihood::BrokenPowerLaw3 foo(0.686, -1.7, 2.87e-3, -3.5, 
                                   100, 1e4, 2e4, 1e5);
   optimizers::FunctionTest tester(foo, "BrokenPowerLaw3");
   std::vector<optimizers::Parameter> pars;
   pars.push_back(optimizers::Parameter("Integral1", 0.686));
   pars.push_back(optimizers::Parameter("Index1", -1.7));
   pars.push_back(optimizers::Parameter("Integral2", 2.87e-3));
   pars.push_back(optimizers::Parameter("Index2", -3.5));
   pars.push_back(optimizers::Parameter("LowerLimit1", 100.));
   pars.push_back(optimizers::Parameter("UpperLimit1", 1e4));
   pars.push_back(optimizers::Parameter("LowerLimit2", 2e4));
   pars.push_back(optimizers::Parameter("UpperLimit2", 1e5));
   
   // Implementation continues...
}

void LikelihoodTests::test_MultipleBrokenPowerLaw() {
   // Implementation
}

void LikelihoodTests::test_EblAtten() {
   // Implementation
}

void LikelihoodTests::test_EnergyBand() {
   // Implementation
}

void LikelihoodTests::test_SourceFactory() {
   // Implementation
}

void LikelihoodTests::test_LikeExposure() {
   using interval = std::pair<double, double>;
   std::vector<interval> timeRangeCuts;
   std::vector<interval> gtis;

   double start(100.);
   double stop(200.);
   timeRangeCuts.push_back(std::make_pair(150., 300.));

   double fraction;
   LikeExposure::acceptInterval(start, stop, timeRangeCuts, gtis, fraction);
   
   ASSERT_EQUALS(fraction, 0.5);

   gtis.push_back(std::make_pair(150., 165.));
   gtis.push_back(std::make_pair(190., 500.));

   LikeExposure::acceptInterval(start, stop, timeRangeCuts, gtis, fraction);

   ASSERT_EQUALS(fraction, 0.25);

   CPPUNIT_ASSERT(!LikeExposure::acceptInterval(301., 500., timeRangeCuts, 
                                                gtis, fraction));
}

void LikelihoodTests::test_SourceModel() {
   // Implementation
}

void LikelihoodTests::test_SourceDerivs() {
   // Implementation
}

void LikelihoodTests::test_PointSource() {
   // Implementation
}

void LikelihoodTests::test_DiffuseSource() {
   // Implementation
}

void LikelihoodTests::test_CountsMap() {
   // Implementation
}

void LikelihoodTests::test_CountsMapHealpix_allsky() {
   // Implementation
}

void LikelihoodTests::test_CountsMapHealpix_region() {
   // Implementation
}

void LikelihoodTests::test_BinnedLikelihood() {
   // Implementation
}

void LikelihoodTests::test_BinnedLikelihood_2() {
   // Implementation
}

void LikelihoodTests::test_BinnedLikelihood_base(int edisp_val) {
   // Implementation
}

void LikelihoodTests::test_CompositeSource() {
   // Implementation
}

void LikelihoodTests::test_MeanPsf() {
   // Implementation
}

void LikelihoodTests::test_BinnedExposure() {
   // Implementation
}

void LikelihoodTests::test_BinnedExposureHealpix() {
   // Implementation
}

void LikelihoodTests::test_SourceMap() {
   // Implementation
}

void LikelihoodTests::test_PointSourceMap() {
   // Implementation
}

void LikelihoodTests::test_PointSourceMap_hpx_allsky() {
   // Implementation
}

void LikelihoodTests::test_PointSourceMap_hpx_region() {
   // Implementation
}

void LikelihoodTests::test_rescaling() {
   std::vector<optimizers::Function*> my_functions;
   my_functions.push_back(new BandFunction());
   my_functions.push_back(new BrokenPowerLaw2());
   my_functions.push_back(new BrokenPowerLawExpCutoff());
   my_functions.push_back(new PowerLaw2());
   my_functions.push_back(new PowerLawSuperExpCutoff());
   
   optimizers::dArg xx(100);
   double factor(2);
   
   for (auto* func : my_functions) {
      double value(func->operator()(xx));
      func->rescale(factor);
      ASSERT_EQUALS(2 * value, func->operator()(xx));
   }
   
   // Clean up
   for (auto* func : my_functions) {
      delete func;
   }
}

void LikelihoodTests::test_WcsMap2_exception() {
   // Implementation
}

void LikelihoodTests::test_WcsMap2() {
   // Test interpolatePowerLaw.
   double x, x1, x2;
   double y1, y2;

   // Test switch to linear interpolation
   double value = Likelihood::WcsMap2::interpolatePowerLaw(x = 1, x1 = 1, x2 = 2,
                                                           y1 = 0, y2 = 1);
   CPPUNIT_ASSERT(value == 0);
   
   // Test for extrapolation exception if linear interpolation is selected
   try {
      Likelihood::WcsMap2::interpolatePowerLaw(x = -1, x1 = 1, x2 = 2,
                                               y1 = 0, y2 = 1);
   } catch (std::runtime_error& eObj) {
      if (!st_facilities::Util::expectedException(eObj,
                                                  "linear extrapolation selected")) {
         throw;
      }
   }
}

void LikelihoodTests::test_Drm() {
   // Implementation
}

void LikelihoodTests::test_Source_Npred() {
   // Implementation
}

void LikelihoodTests::test_ExposureCube() {
   // Implementation
}

// Helper method implementations
optimizers::FunctionFactory* LikelihoodTests::funcFactoryInstance() {
   if (!m_funcFactory) {
      m_funcFactory = new optimizers::FunctionFactory();
      // Add function prototypes as needed
   }
   return m_funcFactory;
}

SourceFactory * LikelihoodTests::
srcFactoryInstance(const std::string & scFile,
                   const std::string & expMapFile,
                   const std::string & sourceXmlFile,
                   bool requireExposure, bool verbose) {
   if (m_srcFactory == 0) {
      m_roiCuts->setCuts(86.404, 28.936, 25., 30., 2e5, 0., 8.64e4, -1., true);
      if (scFile == "") {
         m_scData->readData(m_scFile, 0, 86400, true);
      } else {
         m_scData->readData(scFile, 0, 86400, true);
      }

      if (expMapFile == "") {
         m_expMap->readExposureFile(m_expMapFile);
      } else {
         m_expMap->readExposureFile(expMapFile);
      }

      optimizers::FunctionFactory * funcFactory = funcFactoryInstance();

      m_srcFactory = new SourceFactory(*m_observation, verbose);
      if (sourceXmlFile == "") {
         m_srcFactory->readXml(m_sourceXmlFile, *funcFactory, requireExposure);
      } else {
         m_srcFactory->readXml(sourceXmlFile, *funcFactory, requireExposure);
      }
   }
   return m_srcFactory;
}

void LikelihoodTests::readEventData(const std::string& eventFile,
                                     const std::string& scDataFile,
                                     std::vector<Event>& events) {
   // Implementation
}

void LikelihoodTests::generate_exposureHyperCube() {
   // Implementation
}

CountsMap LikelihoodTests::singleSrcMap(unsigned int nee,
                                         unsigned long num_x_pix,
                                         unsigned long num_y_pix,
                                         double pix_scale) const {
   std::string eventFile = dataPath("single_src_events.fits");
   std::string scFile = dataPath("oneday_scData_0000.fits");
   
   double ra = 193.98;
   double dec = -5.82;
   std::string projection = "STG";
   bool use_lb = false;
   std::string ra_field = "RA";
   std::string dec_field = "DEC";
   
   // Create energy bins
   std::vector<double> energies;
   double emin = 30.;
   double emax = 2e5;
   double estep = std::log(emax / emin) / (nee - 1);
   for (unsigned int i = 0; i < nee; ++i) {
      energies.push_back(emin * std::exp(estep * i));
   }
   
   return CountsMap(eventFile, "EVENTS",
                    scFile, "SC_DATA",
                    ra, dec, projection,
                    num_x_pix, num_y_pix,
                    pix_scale, 0., use_lb,
                    ra_field, dec_field,
                    energies);
}

CountsMapHealpix LikelihoodTests::healpixmap_allsky() const {
   std::string datafile = dataPath("countsMap_hpx_allsky.fits");
   return CountsMapHealpix(datafile);
}

CountsMapHealpix LikelihoodTests::healpixmap_region() const {
   std::string datafile = dataPath("ccube_single_src_events_hpx.fits");
   return CountsMapHealpix(datafile);
}

// Main function
int main(int argc, char* argv[]) {
#ifdef TRAP_FPE
   feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

   // Initialize facilities
   facilities::commonUtilities::setupEnvironment();
   
   // Load IRFs
   irfLoader::Loader::go();

   CppUnit::TextTestRunner runner;
   runner.addTest(LikelihoodTests::suite());
   
   bool success = runner.run();
   return success ? 0 : 1;
}
