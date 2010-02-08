/**
 * @file test.cxx
 * @brief Test program for Likelihood.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/test/test.cxx,v 1.94 2010/02/02 06:48:47 jchiang Exp $
 */

#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <cmath>
#include <cstdio>

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include <cppunit/ui/text/TextTestRunner.h>
#include <cppunit/extensions/HelperMacros.h>

#include "facilities/Util.h"
#include "facilities/commonUtilities.h"

#include "st_facilities/Util.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "optimizers/dArg.h"
#include "optimizers/FunctionFactory.h"
#include "optimizers/FunctionTest.h"
#include "optimizers/Minuit.h"

#include "irfInterface/IrfsFactory.h"
#include "irfInterface/AcceptanceCone.h"
#include "dc1aResponse/loadIrfs.h"

#include "Likelihood/BinnedExposure.h"
#include "Likelihood/BinnedLikelihood.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/DiffRespNames.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/Event.h"
#include "Likelihood/EventContainer.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/FluxBuilder.h"
#include "Likelihood/LikeExposure.h"
#include "Likelihood/LogNormal.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/Observation.h"
#include "Likelihood/PointSource.h"
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

#include "Likelihood/BandFunction.h"
#include "Likelihood/BrokenPowerLaw2.h"
#include "Likelihood/BrokenPowerLawExpCutoff.h"
#include "Likelihood/EblAtten.h"
#include "Likelihood/PowerLaw2.h"
#include "Likelihood/PowerLawSuperExpCutoff.h"
#include "Likelihood/SmoothBrokenPowerLaw.h"

#include "SourceData.h"
#include "XmlDiff.h"

using namespace Likelihood;
using optimizers::Parameter;

class LikelihoodTests : public CppUnit::TestFixture {

   CPPUNIT_TEST_SUITE(LikelihoodTests);

   CPPUNIT_TEST(test_LogNormal);
   CPPUNIT_TEST(test_BandFunction);
   CPPUNIT_TEST(test_SmoothBrokenPowerLaw);
   CPPUNIT_TEST(test_EblAtten);
   CPPUNIT_TEST(test_RoiCuts);
   CPPUNIT_TEST(test_SourceFactory);
   CPPUNIT_TEST(test_XmlBuilders);
   CPPUNIT_TEST(test_LikeExposure);
   CPPUNIT_TEST(test_SourceModel);
   CPPUNIT_TEST(test_SourceDerivs);
   CPPUNIT_TEST(test_PointSource);
   CPPUNIT_TEST(test_DiffuseSource);
   CPPUNIT_TEST(test_CountsMap);
   CPPUNIT_TEST(test_BinnedLikelihood);
   CPPUNIT_TEST(test_MeanPsf);
   CPPUNIT_TEST(test_BinnedExposure);
   CPPUNIT_TEST(test_SourceMap);
   CPPUNIT_TEST(test_rescaling);
   CPPUNIT_TEST(test_DiffRespNames);

   CPPUNIT_TEST_SUITE_END();

public:

   void setUp();
   void tearDown();

   void test_LogNormal();
   void test_BandFunction();
   void test_SmoothBrokenPowerLaw();
   void test_EblAtten();
   void test_RoiCuts();
   void test_SourceFactory();
   void test_XmlBuilders();
   void test_LikeExposure();
   void test_SourceModel();
   void test_SourceDerivs();
   void test_PointSource();
   void test_DiffuseSource();
   void test_CountsMap();
   void test_BinnedLikelihood();
   void test_MeanPsf();
   void test_BinnedExposure();
   void test_SourceMap();
   void test_rescaling();
   void test_DiffRespNames();

private:

   Observation * m_observation;
   static RoiCuts * m_roiCuts;
   static ScData * m_scData;
   static ExposureCube * m_expCube;
   static ExposureMap * m_expMap;
   static ResponseFunctions * m_respFuncs;
   static EventContainer * m_eventCont;

   std::string m_rootPath;
   double m_fracTol;

// File names for m_srcFactory.
   std::string m_scFile;
   std::string m_expMapFile;
   std::string m_sourceXmlFile;

// Test data for m_srcFactory.
   SourceData m_srcData;

// Filenames for test_XmlBuilders.
   std::string m_fluxXmlFile;
   std::string m_srcModelXmlFile;

   optimizers::FunctionFactory * m_funcFactory;
   optimizers::FunctionFactory * funcFactoryInstance();

   SourceFactory * m_srcFactory;
   SourceFactory * srcFactoryInstance(const std::string & scFile="",
                                      const std::string & expMapFile="",
                                      const std::string & sourceXmlFile="",
                                      bool requireExposure=true,
                                      bool verbose=true);

   void readEventData(const std::string & eventFile,
                      const std::string & scDataFile,
                      std::vector<Event> & events);

   void generate_exposureHyperCube();

   CountsMap singleSrcMap(unsigned int nee) const;

   void deleteExpMap();
};

#define ASSERT_EQUALS(X, Y) CPPUNIT_ASSERT(fabs( (X - Y)/Y ) < m_fracTol)

RoiCuts * LikelihoodTests::m_roiCuts(0);
ScData * LikelihoodTests::m_scData(0);
ExposureCube * LikelihoodTests::m_expCube(0);
ExposureMap * LikelihoodTests::m_expMap(0);
ResponseFunctions * LikelihoodTests::m_respFuncs(0);
EventContainer * LikelihoodTests::m_eventCont(0);

void LikelihoodTests::setUp() {
   facilities::commonUtilities::setupEnvironment();
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
   const char * root = facilities::commonUtilities::getPackagePath("Likelihood").c_str();
   if (!root) {  //use relative path from cmt directory
      m_rootPath = "..";
   } else {
      m_rootPath = facilities::commonUtilities::getDataPath("Likelihood");
   }

// Prepare the ResponseFunctions object.
   dc1aResponse::load_irfs();
   irfInterface::IrfsFactory * myFactory 
      = irfInterface::IrfsFactory::instance();
   m_respFuncs->addRespPtr(0, myFactory->create("DC1A::Front"));
   m_respFuncs->addRespPtr(1, myFactory->create("DC1A::Back"));

// Fractional tolerance for double comparisons.
   m_fracTol = 1e-4;

// Use lazy evaluation for m_funcFactory and m_srcFactory.
   m_funcFactory = 0;
   m_srcFactory = 0;

   m_scFile = facilities::commonUtilities::joinPath(m_rootPath,"oneday_scData_0000.fits");
   m_expMapFile = facilities::commonUtilities::joinPath(m_rootPath,"anticenter_expMap.fits");
   m_sourceXmlFile = facilities::commonUtilities::joinPath(m_rootPath, "anticenter_model.xml");
}

void LikelihoodTests::deleteExpMap() {
   delete m_expMap;
   m_expMap = 0;
}

void LikelihoodTests::tearDown() {
   delete m_funcFactory;
   m_funcFactory = 0;
   delete m_srcFactory;
   m_srcFactory = 0;

   delete m_observation;
   m_observation = 0;
 
// @todo Use iterators to traverse RespPtr map for key deletion.
   m_respFuncs->deleteRespPtr(0);
   m_respFuncs->deleteRespPtr(1);
   std::remove(m_fluxXmlFile.c_str());
   std::remove(m_srcModelXmlFile.c_str());
}

void LikelihoodTests::test_LogNormal() {
   Likelihood::LogNormal func(1, 3, 2);
   optimizers::FunctionTest tester(func, "LogNormal");
   std::vector<optimizers::Parameter> params;
   params.push_back(optimizers::Parameter("Prefactor", 1));
   params.push_back(optimizers::Parameter("Log10_Mean", 3));
   params.push_back(optimizers::Parameter("Log10_Sigma", 2));

   std::vector<optimizers::Arg *> args;
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

   std::vector<optimizers::Arg *> args;
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
}

void LikelihoodTests::test_SmoothBrokenPowerLaw() {
   Likelihood::SmoothBrokenPowerLaw foo(10, -2.1, 100, -2.1, 1e3, 0.2);
   optimizers::FunctionTest tester(foo, "SmoothBrokenPowerLaw");
   std::vector<optimizers::Parameter> params;
   params.push_back(optimizers::Parameter("Prefactor", 10));
   params.push_back(optimizers::Parameter("Index1", -2.1));
   params.push_back(optimizers::Parameter("Scale", 100.));
   params.push_back(optimizers::Parameter("Index2", -2.1));
   params.push_back(optimizers::Parameter("BreakValue", 1e3));
   params.push_back(optimizers::Parameter("Beta", 0.2));

   std::vector<optimizers::Arg *> args;
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
}

void LikelihoodTests::test_EblAtten() {
   Likelihood::EblAtten foo(PowerLaw2(1, -2.1, 20, 2e5), 1, 0.5, 0);
   optimizers::FunctionTest tester(foo, "EblAtten::PowerLaw2");
   std::vector<optimizers::Parameter> params;
   params.push_back(optimizers::Parameter("Integral", 1));
   params.push_back(optimizers::Parameter("Index", -2.1));
   params.push_back(optimizers::Parameter("LowerLimit", 20.));
   params.push_back(optimizers::Parameter("UpperLimit", 2e5));
   params.push_back(optimizers::Parameter("tau_norm", 1));
   params.push_back(optimizers::Parameter("redshift", 0.5));
   params.push_back(optimizers::Parameter("ebl_model", 0));

   std::vector<optimizers::Arg *> args;
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
}

void LikelihoodTests::test_RoiCuts() {
   m_roiCuts->setCuts(193.98, -5.82, 20, 30, 3.1623e5, 0, 1e12, -1, true);

// Compare to known default values.
   std::vector<std::pair<double, double> > tlims;
   m_roiCuts->getTimeCuts(tlims);
   static double tmin = 0;
   static double tmax = 1e12;
   CPPUNIT_ASSERT(fabs(tlims[0].first - tmin) == 0);
   ASSERT_EQUALS(tlims[0].second, tmax);

   std::pair<double, double> energies = m_roiCuts->getEnergyCuts();
   static double emin = 30.;
   static double emax = 3.1623e5;
   ASSERT_EQUALS(energies.first, emin);
   ASSERT_EQUALS(energies.second, emax);

   static double ra = 193.98;
   static double dec = -5.82;
   static double radius = 20.;
   irfInterface::AcceptanceCone roiCone(astro::SkyDir(ra, dec), radius);
   CPPUNIT_ASSERT(roiCone == m_roiCuts->extractionRegion());
   double my_ra, my_dec;
   m_roiCuts->getRaDec(my_ra, my_dec);
   ASSERT_EQUALS(my_ra, ra);
   ASSERT_EQUALS(my_dec, dec);
}

void LikelihoodTests::test_SourceFactory() {

   SourceFactory * srcFactory = srcFactoryInstance();

   std::vector<std::string> srcNames;
   m_srcData.getSrcNames(srcNames);

   for (int i = 0; i < 5; i++) {
      std::string name = srcNames[i];
      Source * src = srcFactory->create(name);
      Source::FuncMap srcFuncs = src->getSrcFuncs();
      CPPUNIT_ASSERT(src->getType() == m_srcData.srcType(name));
      CPPUNIT_ASSERT(srcFuncs.count("Spectrum"));
      const std::vector<std::string> & paramNames = m_srcData.paramNames();
      for (unsigned int j = 0; j < paramNames.size(); j++) {
         ASSERT_EQUALS(srcFuncs["Spectrum"]->getParamValue(paramNames[j]),
                       m_srcData.paramObject(name, paramNames[j]).getValue());
      }
      if (m_srcData.srcType(name) == "Diffuse") {
         CPPUNIT_ASSERT(srcFuncs.count("SpatialDist"));
         CPPUNIT_ASSERT(srcFuncs["SpatialDist"]->genericName() 
                == m_srcData.spatialModel(name));
      } else {
         CPPUNIT_ASSERT(srcFuncs.count("Position"));
         CPPUNIT_ASSERT(srcFuncs["Position"]->genericName() 
                == m_srcData.spatialModel(name));
      }
   }
}

void LikelihoodTests::test_XmlBuilders() {

   SourceFactory * srcFactory = srcFactoryInstance();
   std::vector<std::string> srcNames;
   srcFactory->fetchSrcNames(srcNames);

   FluxBuilder fluxBuilder(30, 2e5);
   SourceModelBuilder srcModelBuilder("", "source library");
   for (unsigned int i = 0; i < srcNames.size(); i++) {
      Source * src = srcFactory->create(srcNames[i]);
      fluxBuilder.addSource(*src);
      srcModelBuilder.addSource(*src);
   }
   m_fluxXmlFile = facilities::commonUtilities::joinPath(m_rootPath, "fluxBuilder.xml");
   fluxBuilder.write(m_fluxXmlFile);

   m_srcModelXmlFile = facilities::commonUtilities::joinPath(m_rootPath, "srcModelBuilder.xml");
   srcModelBuilder.write(m_srcModelXmlFile);

   XmlDiff xmlDiff(m_sourceXmlFile, m_srcModelXmlFile, "source", "name");
   CPPUNIT_ASSERT(xmlDiff.compare());
}

void LikelihoodTests::test_LikeExposure() {

   typedef std::pair<double, double> interval;
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

   SourceFactory * srcFactory = srcFactoryInstance();
   std::vector<std::string> srcNames;
   srcFactory->fetchSrcNames(srcNames);

   SourceModel srcModel(*m_observation);

   for (unsigned int i = 0; i < srcNames.size(); i++) {
      srcModel.addSource(srcFactory->create(srcNames[i]));
   }

// Test the parameter values contained in srcModel.
   m_srcData.getSrcNames(srcNames);
   const std::vector<std::string> & paramNames = m_srcData.paramNames();
   for (unsigned int i = 0; i < srcNames.size(); i++) {
      std::string name = srcNames[i];
      for (unsigned int j = 0; j < 3; j++) {
         ASSERT_EQUALS(srcModel.getParamValue(paramNames[j], "Spectrum", name),
                       m_srcData.paramObject(name, paramNames[j]).getValue());
      }
   }

// Check free Parameters.
   srcModel.getSrcNames(srcNames);
   std::vector<double> my_parValues;
   for (unsigned int i = 0; i < srcNames.size(); i++) {
      std::string name = srcNames[i];
      for (unsigned int j = 0; j < 3; j++) {
         if (m_srcData.paramObject(name, paramNames[j]).isFree()) {
            my_parValues.push_back(
               m_srcData.paramObject(name, paramNames[j]).getValue());
         }
      }
   }
   std::vector<double> freeParValues;
   srcModel.getFreeParamValues(freeParValues);
   CPPUNIT_ASSERT(my_parValues.size() == freeParValues.size());
   for (unsigned int i = 0; i < my_parValues.size(); i++) {
      ASSERT_EQUALS(my_parValues[i], freeParValues[i]);
   }

   std::vector<double> saved_values;
   saved_values = freeParValues;
   std::vector<double> newFreeParams;
   for (unsigned int i = 0; i < freeParValues.size(); i++) {
      newFreeParams.push_back(freeParValues[i]*1.1);
   }
   srcModel.setFreeParamValues(newFreeParams);
   srcModel.getFreeParamValues(freeParValues);
   CPPUNIT_ASSERT(newFreeParams.size() == freeParValues.size());
   for (unsigned int i = 0; i < newFreeParams.size(); i++) {
      ASSERT_EQUALS(newFreeParams[i], freeParValues[i]);
   }

// Reset the values.
   srcModel.setFreeParamValues(saved_values);

// Derivative tests.  Although the SourceModel::value(Arg &) and
// SourceModel::derivByParam(Arg &, ...) methods are largely
// meaningless, those interfaces provide a way of checking the
// delegation to the corresponding methods in the Function objects
// that make-up the individual Source components.

// Note that we have to delete the Galactic Diffuse source since the
// SpatialMap class, used as a functor by the Galactic Diffuse model,
// will try to downcast the dArg variable to a SkyDirArg.
   srcModel.deleteSource("Galactic Diffuse");

// Delete Extragalactic Diffuse since its ConstantValue function has a
// much larger value than the other functions that the Source model
// value(...) method sums together, unity vs O(1e-10), making the
// numerical derivatives difficult to compute.
   srcModel.deleteSource("Extragalactic Diffuse");

// Consider the derivatives wrt free parameters:
   optimizers::dArg x(1000.);
   std::vector<double> params_save;
   srcModel.getFreeParamValues(params_save);
   std::vector<double> sm_derivs;
   srcModel.getFreeDerivs(x, sm_derivs);
   std::vector<double> params = params_save;

   for (unsigned int i = 0; i < sm_derivs.size(); i++) {
// Compute the numerical derivative wrt this parameter.
      double delta_param = fabs(params_save[i]/1e7);
      double num_deriv = srcModel(x);
      params[i] += delta_param;
      srcModel.setFreeParamValues(params);
      num_deriv -= srcModel(x);
      num_deriv /= delta_param;
      ASSERT_EQUALS(sm_derivs[i], -num_deriv);

// Reset the parameters for next time around.
      srcModel.setFreeParamValues(params_save);
      params = params_save;
   }

// Derivatives wrt all parameters.
   srcModel.getParamValues(params_save);
   srcModel.getDerivs(x, sm_derivs);
   params = params_save;

   for (unsigned int i = 0; i < sm_derivs.size(); i++) {
// Compute the numerical derivative wrt this parameter.
      double delta_param = fabs(params_save[i]/1e7);
      double num_deriv = srcModel(x);
      params[i] += delta_param;
      srcModel.setParamValues(params);
      num_deriv -= srcModel(x);
      num_deriv /= delta_param;
// We cannot use ASSERT_EQUALS here since some of the derivatives
// (e.g., those wrt ra, dec) are identically zero.
      CPPUNIT_ASSERT(fabs(sm_derivs[i] + num_deriv) < m_fracTol);

// Reset the parameters for next time around.
      srcModel.setParamValues(params_save);
      params = params_save;
   }
}

void LikelihoodTests::test_SourceDerivs() {
   SourceFactory * srcFactory = srcFactoryInstance();
   CPPUNIT_ASSERT(srcFactory != 0);
   std::vector<std::string> srcNames;
   srcFactory->fetchSrcNames(srcNames);

// Create a 1 GeV Front-converting event from the anticenter at
// MET=1000 s.
   double ra(86.4);
   double dec(28.9);
   double energy(1e3);
   double time(1e3);
   CPPUNIT_ASSERT(m_scData != 0);
   astro::SkyDir zAxis = m_scData->zAxis(time);
   double muZenith(1.);
   int type = 0;
   Event myEvent(ra, dec, energy, time, zAxis.ra(), zAxis.dec(), 
                 muZenith, m_respFuncs->useEdisp(), m_respFuncs->respName(), 
                 type);

// Loop over Sources and check fluxDensity and Npred derivatives
// for each.
//   for (int i = 0; i < 3; i++) {
   for (int i = 0; i < 2; i++) {
      Source * src = srcFactory->create(srcNames[i]);
      CPPUNIT_ASSERT(src != 0);
      if (src->getType() == "Diffuse") {
         DiffuseSource *diffuse_src = dynamic_cast<DiffuseSource *>(src);
         myEvent.computeResponse(*diffuse_src, *m_respFuncs);
      }
      Source::FuncMap srcFuncs = src->getSrcFuncs();
      CPPUNIT_ASSERT(srcFuncs.count("Spectrum"));

      std::vector<double> paramValues;
      srcFuncs["Spectrum"]->getFreeParamValues(paramValues);

      std::vector<std::string> paramNames;
      srcFuncs["Spectrum"]->getFreeParamNames(paramNames);
      CPPUNIT_ASSERT(paramValues.size() == paramNames.size());

      double fluxDensity0 = src->fluxDensity(myEvent);
      double Npred0 = src->Npred();
      double eps = 1e-7;
      for (unsigned int j = 0; j < paramNames.size(); j++) {
         double delta = eps*paramValues[j];
         srcFuncs["Spectrum"]->setParam(paramNames[j], paramValues[j]+delta);

         double fluxDensity1 = src->fluxDensity(myEvent);
         double fluxDensityDeriv_est = (fluxDensity1 - fluxDensity0)/delta;
         double fluxDensityDeriv 
            = src->fluxDensityDeriv(myEvent, paramNames[j]);
         ASSERT_EQUALS(fluxDensityDeriv_est, fluxDensityDeriv);

         double Npred1 = src->Npred();
         double NpredDeriv_est = (Npred1 - Npred0)/delta;
         double NpredDeriv = src->NpredDeriv(paramNames[j]);
         ASSERT_EQUALS(NpredDeriv_est, NpredDeriv);

         srcFuncs["Spectrum"]->setParam(paramNames[j], paramValues[j]);
      }
   }
}

void LikelihoodTests::test_PointSource() {
   std::string eventFile = facilities::commonUtilities::joinPath(m_rootPath, "single_src_events_0000.fits");

   tearDown();
   setUp();

   std::vector<Event> events;
   readEventData(eventFile, m_scFile, events);

   SourceFactory * srcFactory = srcFactoryInstance();

   Source * src = srcFactory->create("Crab Pulsar");

   double Nobs = 0;
   for (unsigned int i = 0; i < events.size(); i++) {
      if (m_roiCuts->accept(events[i])) {
         Nobs++;
      }
   }
   dynamic_cast<PointSource *>(src)->setDir(83.57, 22.01, true, false);
   std::ostringstream debug_output;
   debug_output << Nobs << "  "
                << src->Npred() << "\n";

// Consider the observation over 0.1 day intervals over the one day,
// resetting the ROI accordingly and force the PointSource exposure to
// be recomputed for each interval.
//
// Also test that resetting the RoiCuts energy bounds results in
// a proper PointSource exposure recalculation.
   double eminVals[] = {30., 40., 50., 60., 70.,
                        80., 100., 150., 200., 300.};
   double chi2 = 0;
   double tstep = 8.64e4/10.;
   for (int j = 0; j < 10; j++) {
      double tmin = j*tstep;
      double tmax = tmin + tstep;
      m_roiCuts->setCuts(86.4, 28.9, 25., eminVals[j], 2e5, tmin, tmax, -1.,
                         true);

      dynamic_cast<PointSource *>(src)->setDir(83.57, 22.01, true, false);

      Nobs = 0;
      for (unsigned int i = 0; i < events.size(); i++) {
         if (m_roiCuts->accept(events[i])) {
            Nobs++;
         }
      }
      double Npred = src->Npred();
      chi2 += pow((Nobs - Npred), 2)/Nobs;
      debug_output << j << "  " 
                   << Nobs << "  "
                   << Npred << "\n";
   }
   debug_output << "chi^2 = " << chi2;
   if (chi2 >= 20.) {
      std::cout << debug_output.str() << std::endl;
   }
   CPPUNIT_ASSERT(chi2 < 20.);
}

void LikelihoodTests::test_DiffuseSource() {
   std::string eventFile = facilities::commonUtilities::joinPath(m_rootPath, "galdiffuse_events_0000.fits");

   std::vector<Event> events;
   readEventData(eventFile, m_scFile, events);

   astro::SkyDir anticenter(180., 0., astro::SkyDir::GALACTIC);

   std::ostringstream debug_output;

   double chi2 = 0;
   for (int i = 0; i < 5; i++) {
      tearDown();
      deleteExpMap();
      setUp();

      std::ostringstream expMapFile;
      expMapFile << m_rootPath << "/expMap_" << i << ".fits";

      SourceFactory * srcFactory = srcFactoryInstance("", expMapFile.str());

      double tstep = 0.2*8.64e4;
      double tmin = i*tstep;
      double tmax = tmin + tstep;
      m_roiCuts->setCuts(anticenter.ra(), anticenter.dec(), 20., 30., 2e5,
                         tmin, tmax, -1., true);

      Source * src = srcFactory->create("Galactic Diffuse");

      double Nobs = 0;
      for (unsigned int j = 0; j < events.size(); j++) {
         if (m_roiCuts->accept(events[j])) {
            Nobs++;
         }
      }
      double Npred = src->Npred();
      chi2 += pow((Nobs - Npred), 2)/Nobs;
      debug_output << i << "  " 
                << Nobs << "  "
                << Npred << "  " 
                << expMapFile.str() << std::endl;
      delete src;
   }
   debug_output << "chi^2 = " << chi2 << std::endl;
   if (chi2 >= 6.) {
      std::cout << debug_output.str() << std::endl;
   }
   CPPUNIT_ASSERT(chi2 < 6.);
}

void LikelihoodTests::generate_exposureHyperCube() {
   srcFactoryInstance();
   std::vector<std::pair<double, double> > timeCuts;
   m_roiCuts->getTimeCuts(timeCuts);
   std::vector< std::pair<double, double> > gtis;
   m_roiCuts->getGtis(gtis);
   LikeExposure exposure(1., 0.025, timeCuts, timeCuts);
   const tip::Table * scData = tip::IFileSvc::instance().readTable(m_scFile,
                                                                   "SC_DATA");
   exposure.load(scData, false);
   std::string output_file = facilities::commonUtilities::joinPath(m_rootPath, "/expcube_1_day.fits");
   exposure.write(output_file);
}

CountsMap LikelihoodTests::singleSrcMap(unsigned int nee) const {
   std::string eventFile = facilities::commonUtilities::joinPath(m_rootPath, "single_src_events_0000.fits");
   double ra(83.57);
   double dec(22.01);
   unsigned long npts(40);
   double emin(30.);
   double emax(2e5);
   CountsMap dataMap(eventFile, "EVENTS", m_scFile, "SC_DATA", 
                     ra, dec, "CAR", npts, npts,
                     0.25, 0, false, "RA", "DEC", emin, emax, nee);
   const tip::Table * events 
      = tip::IFileSvc::instance().readTable(eventFile, "events");
   dataMap.binInput(events->begin(), events->end());
   return dataMap;
}

void LikelihoodTests::test_CountsMap() {
   CountsMap dataMap(singleSrcMap(21));
   dataMap.writeOutput("test_CountsMap", "countsMap.fits");
   CountsMap dataMap2("countsMap.fits");
   for (unsigned int i = 0; i < dataMap.data().size(); i++) {
      CPPUNIT_ASSERT(dataMap.data()[i] == dataMap2.data()[i]);
   }
}

void LikelihoodTests::test_BinnedLikelihood() {
   std::string exposureCubeFile = facilities::commonUtilities::joinPath(m_rootPath, "expcube_1_day.fits");
   if (!st_facilities::Util::fileExists(exposureCubeFile)) {
      generate_exposureHyperCube();
   }
   m_expCube->readExposureCube(exposureCubeFile);

   SourceFactory * srcFactory = srcFactoryInstance();
   (void)(srcFactory);

   CountsMap dataMap(singleSrcMap(21));

   BinnedLikelihood binnedLogLike(dataMap, *m_observation);
   std::string Crab_model = facilities::commonUtilities::joinPath(m_rootPath, "Crab_model.xml");
   binnedLogLike.readXml(Crab_model, *m_funcFactory);
   binnedLogLike.saveSourceMaps("srcMaps.fits");

// Try to fit using binned model.
   optimizers::Minuit my_optimizer(binnedLogLike);
   int verbose(0);
   double tol(1e-5);
   my_optimizer.find_min(verbose, tol, optimizers::RELATIVE);

    std::vector<double> params;
//    binnedLogLike.getFreeParamValues(params);
//    for (unsigned int i = 0; i < params.size(); i++) {
//       std::cout << "parameter " << i << ": " << params[i] << "\n";
//    }
//    std::cout << std::endl;

   CountsMap * modelMap = binnedLogLike.createCountsMap();

   dataMap.writeOutput("test_Likelihood", "dataMap.fits");
   modelMap->writeOutput("test_Likelihood", "modelMap.fits");

   const std::vector<float> & data = dataMap.data();
   double dataSum(0);
   for (unsigned int i = 0; i < data.size(); i++) {
      dataSum += data[i];
   }
//   std::cout << "Total counts in data map: " << dataSum << std::endl;

   const std::vector<float> & model = modelMap->data();
   double modelSum(0);
   for (unsigned int i = 0; i < model.size(); i++) {
      modelSum += model[i];
   }
//   std::cout << "Total model counts: " << modelSum << std::endl;

   CPPUNIT_ASSERT(fabs(modelSum - dataSum)/dataSum < 1e-2);

   unsigned long npts = dataMap.imageDimension(0);
   std::vector<double> energies;
   modelMap->getAxisVector(2, energies);
   for (unsigned int i = 0; i < energies.size()-1; i++) {
      double data_counts(0);
      double model_counts(0);
      for (unsigned int j = 0; j < npts*npts; j++) {
         int indx = npts*npts*i + j;
         data_counts += data[indx];
         model_counts += model[indx];
      }
//       std::cout << energies[i] << "  "
//                 << data_counts << "  "
//                 << model_counts << "\n";
   }
//    std::cout << std::endl;

   std::vector<double> derivs;
   binnedLogLike.getFreeDerivs(derivs);
   binnedLogLike.getFreeParamValues(params);

//   std::cout << "Testing derivatives" << std::endl;
   double logLike0 = binnedLogLike.value();
   double eps(1e-7);
   for (unsigned int i = 0; i < params.size(); i++) {
      std::vector<double> new_params = params;
      double delta = eps*new_params[i];
      new_params[i] += delta;
      binnedLogLike.setFreeParamValues(new_params);
      double logLike = binnedLogLike.value();
//       std::cout << i << "  ";
//       std::cout << derivs[i] << "  ";
//       std::cout << logLike << "  " << logLike0 << "  ";
//       std::cout << (logLike - logLike0)/delta << std::endl;
      
// Another weak test.
      double num_deriv = fabs((derivs[i] - (logLike - logLike0)/delta)
                              /derivs[i]);
//       std::cout << "numerical deriv: " << num_deriv << std::endl;
      CPPUNIT_ASSERT(num_deriv < 5e-2);
   }

   delete modelMap;
}

void LikelihoodTests::test_MeanPsf() {
   std::string exposureCubeFile = facilities::commonUtilities::joinPath(m_rootPath, "expcube_1_day.fits");
   if (!st_facilities::Util::fileExists(exposureCubeFile)) {
      generate_exposureHyperCube();
   }
   m_expCube->readExposureCube(exposureCubeFile);

   SourceFactory * srcFactory = srcFactoryInstance();
   (void)(srcFactory);

   double ee[] = {1e2, 1e3, 1e4};
   std::vector<double> energies(ee, ee+3);

   MeanPsf Crab_psf(83.57, 22.01, energies, *m_observation);
   int npts(200);
   double tstep = log(70./1e-4)/(npts-1.);
   std::vector<double> thetas;
   for (int i = 0; i < npts; i++) {
      thetas.push_back(1e-2*exp(i*tstep)*M_PI/180.);
   }

   for (unsigned int k = 0; k < energies.size(); k++) {
      std::vector<double> integrand;
      for (int i = 0; i < npts; i++) {
         double psf = Crab_psf(energies[k], thetas[i]*180./M_PI);
         integrand.push_back(psf*sin(thetas[i])*2.*M_PI);
      }
      TrapQuad my_trap(thetas, integrand);
//       std::cout << my_trap.integral() << std::endl;
// Yes, this test is pretty weak.
      CPPUNIT_ASSERT(fabs(my_trap.integral() - 1.) < 0.03);
   }
}

void LikelihoodTests::test_BinnedExposure() {
   std::string exposureCubeFile = facilities::commonUtilities::joinPath(m_rootPath, "expcube_1_day.fits");
   if (!st_facilities::Util::fileExists(exposureCubeFile)) {
      generate_exposureHyperCube();
   }
   m_expCube->readExposureCube(exposureCubeFile);

   SourceFactory * srcFactory = srcFactoryInstance();
   (void)(srcFactory);

   std::vector<double> energies;
   unsigned int npts(20);
   double emin(20.);
   double emax(2e5);
   double estep = log(emax/emin)/(npts-1.);
   for (unsigned int i = 0; i < npts; i++) {
      energies.push_back(emin*exp(i*estep));
   }
   BinnedExposure binnedExposure(energies, *m_observation);

   std::string filename("binnedExposure.fits");

   binnedExposure.writeOutput(filename);

   BinnedExposure map2(filename);

   double ra(180.);
   double dec(0.);
   for (unsigned int i = 0; i < npts; i++) {
      ASSERT_EQUALS(binnedExposure(energies[i], ra, dec),
                    map2(energies[i], ra, dec));
   }
}

void LikelihoodTests::test_SourceMap() {
   std::string exposureCubeFile = facilities::commonUtilities::joinPath(m_rootPath, "expcube_1_day.fits");
   if (!st_facilities::Util::fileExists(exposureCubeFile)) {
      generate_exposureHyperCube();
   }
   m_expCube->readExposureCube(exposureCubeFile);
   
   CountsMap dataMap(singleSrcMap(3));
   dataMap.writeOutput("test.cxx", "cntsMap.fits");
   
   SourceFactory * srcFactory = srcFactoryInstance("", "", "", false);
 //   Source * src = srcFactory->create("Galactic Diffuse");
   Source * src =  srcFactory->create("Crab Pulsar");

   SourceMap srcMap(src, &dataMap, *m_observation);
}

void LikelihoodTests::test_rescaling() {
   std::vector<optimizers::Function *> my_functions;
   my_functions.push_back(new BandFunction());
   my_functions.push_back(new BrokenPowerLaw2());
   my_functions.push_back(new BrokenPowerLawExpCutoff());
   my_functions.push_back(new PowerLaw2());
   my_functions.push_back(new PowerLawSuperExpCutoff());
   optimizers::dArg xx(100);
   double factor(2);
   for (size_t i(0); i < my_functions.size(); i++) {
      double value(my_functions.at(i)->operator()(xx));
      my_functions.at(i)->rescale(factor);
      ASSERT_EQUALS(2*value, my_functions.at(i)->operator()(xx));
   }
}

void LikelihoodTests::test_DiffRespNames() {
   DiffRespNames foo;
   
   foo.addColumn("P6_V1_DIFFUSE__Extragalactic Diffuse");
   foo.addColumn("P6_V1_DIFFUSE__GALPROP Diffuse");
   
   CPPUNIT_ASSERT("P6_V1_DIFFUSE__Extragalactic Diffuse" == foo[0]);
   CPPUNIT_ASSERT("P6_V1_DIFFUSE__Extragalactic Diffuse" == foo["DIFRSP0"]);
                 
   CPPUNIT_ASSERT("P6_V1_DIFFUSE__GALPROP Diffuse" == foo[1]);
   CPPUNIT_ASSERT("P6_V1_DIFFUSE__GALPROP Diffuse" == foo["DIFRSP1"]);

   CPPUNIT_ASSERT("DIFRSP0" == foo.key("P6_V1_DIFFUSE__Extragalactic Diffuse"));
   CPPUNIT_ASSERT("DIFRSP1" == foo.key("P6_V1_DIFFUSE__GALPROP Diffuse"));
   
   try {
      std::cout << foo[2] << std::endl;
   } catch (std::out_of_range & eObj) {
   }
   try {
      std::cout << foo["DIFRSP2"] << std::endl;
   } catch (std::out_of_range & eObj) {
   }
   try {
      std::cout << foo["DIFRPSs0"] << std::endl;
   } catch (DiffRespNameError & eObj) {
   }
}

void LikelihoodTests::readEventData(const std::string &eventFile,
                                    const std::string &scDataFile,
                                    std::vector<Event> &events) {
   events.clear();
   m_scData->readData(scDataFile, 0, 86400, true);

   tip::Table * eventTable = 
      tip::IFileSvc::instance().editTable(eventFile, "events");

   double ra, dec, energy, time(0), zenith_angle;
   int conversion_layer, type;

   tip::Table::Iterator it = eventTable->begin();
   tip::Table::Record & row = *it;
   
   for ( ; it != eventTable->end(); ++it) {
      astro::SkyDir zAxis = m_scData->zAxis(time);
      row["conversion_layer"].get(conversion_layer);
      if (conversion_layer < 12) {
         type = 0;
      } else {
         type = 1;
      }
      row["ra"].get(ra);
      row["dec"].get(dec);
      row["energy"].get(energy);
      row["time"].get(time);
      row["zenith_angle"].get(zenith_angle);
      events.push_back(Event(ra, dec, energy, time, zAxis.ra(), zAxis.dec(),
                             cos(zenith_angle*M_PI/180.),
                             m_respFuncs->useEdisp(),
                             m_respFuncs->respName(), type));
   }
}   

optimizers::FunctionFactory * LikelihoodTests::funcFactoryInstance() {
   if (m_funcFactory == 0) {
      m_funcFactory = new optimizers::FunctionFactory();
      m_funcFactory->addFunc("SkyDirFunction", new SkyDirFunction(), false);
      m_funcFactory->addFunc("SpatialMap", new SpatialMap(), false);
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

int main(int iargc, char * argv[]) {

#ifdef TRAP_FPE
   feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif

   if (iargc > 1 && std::string(argv[1]) == "-d") { // debug mode
      LikelihoodTests testObj;
      testObj.setUp();
      testObj.test_LogNormal();
      testObj.tearDown();

      testObj.setUp();
      testObj.test_BandFunction();
      testObj.tearDown();

      testObj.setUp();
      testObj.test_EblAtten();
      testObj.tearDown();

      testObj.setUp();
      testObj.test_RoiCuts();
      testObj.tearDown();

      testObj.setUp();
      testObj.test_SourceFactory();
      testObj.tearDown();

      testObj.setUp();
      testObj.test_XmlBuilders();
      testObj.tearDown();

      testObj.setUp();
      testObj.test_LikeExposure();
      testObj.tearDown();

      testObj.setUp();
      testObj.test_SourceModel();
      testObj.tearDown();

      testObj.setUp();
      testObj.test_SourceDerivs();
      testObj.tearDown();

      testObj.setUp();
      testObj.test_PointSource();
      testObj.tearDown();

      testObj.setUp();
      testObj.test_DiffuseSource();
      testObj.tearDown();

      testObj.setUp();
      testObj.test_CountsMap();
      testObj.tearDown();

      testObj.setUp();
      testObj.test_BinnedLikelihood();
      testObj.tearDown();

      testObj.setUp();
      testObj.test_MeanPsf();
      testObj.tearDown();

      testObj.setUp();
      testObj.test_BinnedExposure();
      testObj.tearDown();

      testObj.setUp();
      testObj.test_SourceMap();
      testObj.tearDown();

      testObj.setUp();
      testObj.test_rescaling();
      testObj.tearDown();
   } else {
      CppUnit::TextTestRunner runner;
      runner.addTest(LikelihoodTests::suite());
      bool result = runner.run();
      if (!result) return 1;
   }
}
