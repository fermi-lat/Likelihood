/**
 * @file test.cxx
 * @brief Test program for Likelihood.  Use CppUnit-like idioms.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/test/test.cxx,v 1.9 2004/02/29 03:00:29 jchiang Exp $
 */

#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <cmath>
#include <cassert>
#include <cstdio>

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include "facilities/Util.h"

#include "tuple/ITable.h"

#include "optimizers/dArg.h"
#include "optimizers/FunctionFactory.h"

#include "latResponse/AcceptanceCone.h"
#include "latResponse/IrfsFactory.h"

#include "Likelihood/DiffuseSource.h"
#include "Likelihood/Event.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/FluxBuilder.h"
#include "Likelihood/SourceModelBuilder.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/ScData.h"
#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/Source.h"
#include "Likelihood/SourceFactory.h"
#include "Likelihood/SourceModel.h"
#include "Likelihood/SpatialMap.h"

#include "SourceData.h"
#include "XmlDiff.h"

using namespace Likelihood;
using optimizers::Parameter;

class LikelihoodTests {

public:

   LikelihoodTests() {
      setUp();
   }

   ~LikelihoodTests() {
      tearDown();
   }

   void setUp();
   void tearDown();

   void test_RoiCuts();
   void test_SourceFactory();
   void test_XmlBuilders();
   void test_SourceModel();
   void test_SourceDerivs();
   void test_PointSource();
   void test_DiffuseSource();

private:

   std::string m_rootPath;
   double m_fracTol;

// File names for m_srcFactory.
   std::string m_roiFile;
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
   SourceFactory * srcFactoryInstance(const std::string & roiFile="",
                                      const std::string & scFile="",
                                      const std::string & expMapFile="",
                                      const std::string & sourceXmlFile="",
                                      bool verbose=false);

   void readEventData(const std::string & eventFile,
                      const std::string & scDataFile,
                      std::vector<Event> & events);
   
};

#define ASSERT_EQUALS(X, Y) assert(fabs( (X - Y)/Y ) < m_fracTol)

void LikelihoodTests::setUp() {
// Get root path to test data.
   const char * root = ::getenv("LIKELIHOODROOT");
   if (!root) {  //use relative path from cmt directory
      m_rootPath = "..";
   } else {
      m_rootPath = std::string(root);
   }
// Prepare the ResponseFunctions object.
   latResponse::IrfsFactory irfsFactory;
   ResponseFunctions::addRespPtr(2, irfsFactory.create("DC1::Front"));
   ResponseFunctions::addRespPtr(3, irfsFactory.create("DC1::Back"));
//    ResponseFunctions::addRespPtr(2, irfsFactory.create("Glast25::Front"));
//    ResponseFunctions::addRespPtr(3, irfsFactory.create("Glast25::Back"));   
   
// Fractional tolerance for double comparisons.
   m_fracTol = 1e-4;

// Use lazy evaluation for m_funcFactory and m_srcFactory.
   m_funcFactory = 0;
   m_srcFactory = 0;

   m_roiFile = m_rootPath + "/data/anticenter_Roi.xml";
   m_scFile = m_rootPath + "/data/oneday_scData_0000.fits";
   m_expMapFile = m_rootPath + "/data/anticenter_expMap.fits";
   m_sourceXmlFile = m_rootPath + "/data/anticenter_model.xml";
}

void LikelihoodTests::tearDown() {
   delete m_funcFactory;
   delete m_srcFactory;
// @todo Use iterators to traverse RespPtr map for key deletion.
   ResponseFunctions::deleteRespPtr(2);
   ResponseFunctions::deleteRespPtr(3);
   std::remove(m_fluxXmlFile.c_str());
   std::remove(m_srcModelXmlFile.c_str());
}

void LikelihoodTests::test_RoiCuts() {
   RoiCuts * roiCuts = RoiCuts::instance();
   roiCuts->setCuts();
   std::string xmlFile = m_rootPath + "/data/myROI.xml";
   roiCuts->writeXml(xmlFile, "my_ROI");
   RoiCuts::setCuts(xmlFile);

// Compare to known default values.
   std::vector<std::pair<double, double> > tlims;
   roiCuts->getTimeCuts(tlims);
   static double tmin = 0;
   static double tmax = 1e12;
   assert(fabs(tlims[0].first - tmin) == 0);
   ASSERT_EQUALS(tlims[0].second, tmax);

   std::pair<double, double> energies = roiCuts->getEnergyCuts();
   static double emin = 30.;
   static double emax = 3.1623e5;
   ASSERT_EQUALS(energies.first, emin);
   ASSERT_EQUALS(energies.second, emax);

   static double ra = 193.98;
   static double dec = -5.82;
   static double radius = 20.;
   latResponse::AcceptanceCone roiCone(astro::SkyDir(ra, dec), radius);
   assert(roiCone == roiCuts->extractionRegion());
   double my_ra, my_dec;
   RoiCuts::getRaDec(my_ra, my_dec);
   ASSERT_EQUALS(my_ra, ra);
   ASSERT_EQUALS(my_dec, dec);

   std::remove(xmlFile.c_str());
   std::cout << ".";
}

void LikelihoodTests::test_SourceFactory() {

   SourceFactory * srcFactory = srcFactoryInstance();

   std::vector<std::string> srcNames;
   m_srcData.getSrcNames(srcNames);

   for (int i = 0; i < 5; i++) {
      std::string name = srcNames[i];
      Source * src = srcFactory->create(name);
      Source::FuncMap srcFuncs = src->getSrcFuncs();
      assert(src->getType() == m_srcData.srcType(name));
      assert(srcFuncs.count("Spectrum"));
      const std::vector<std::string> & paramNames = m_srcData.paramNames();
      for (unsigned int j = 0; j < paramNames.size(); j++) {
         ASSERT_EQUALS(srcFuncs["Spectrum"]->getParamValue(paramNames[j]),
                       m_srcData.paramObject(name, paramNames[j]).getValue());
      }
      if (m_srcData.srcType(name) == "Diffuse") {
         assert(srcFuncs.count("SpatialDist"));
         assert(srcFuncs["SpatialDist"]->genericName() 
                == m_srcData.spatialModel(name));
      } else {
         assert(srcFuncs.count("Position"));
         assert(srcFuncs["Position"]->genericName() 
                == m_srcData.spatialModel(name));
      }
   }
   std::cout << ".";
}

void LikelihoodTests::test_XmlBuilders() {

   SourceFactory * srcFactory = srcFactoryInstance();
   std::vector<std::string> srcNames;
   srcFactory->fetchSrcNames(srcNames);

   FluxBuilder fluxBuilder;
   SourceModelBuilder srcModelBuilder("", "source library");
   for (unsigned int i = 0; i < srcNames.size(); i++) {
      Source * src = srcFactory->create(srcNames[i]);
      fluxBuilder.addSource(*src);
      srcModelBuilder.addSource(*src);
   }
   m_fluxXmlFile = m_rootPath + "/data/fluxBuilder.xml";
   fluxBuilder.write(m_fluxXmlFile);

   m_srcModelXmlFile = m_rootPath + "/data/srcModelBuilder.xml";
   srcModelBuilder.write(m_srcModelXmlFile);

   XmlDiff xmlDiff(m_sourceXmlFile, m_srcModelXmlFile, "source", "name");
   assert(xmlDiff.compare());

   std::cout << ".";
}

void LikelihoodTests::test_SourceModel() {

   SourceFactory * srcFactory = srcFactoryInstance();
   std::vector<std::string> srcNames;
   srcFactory->fetchSrcNames(srcNames);
   
   SourceModel srcModel;

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
   assert(my_parValues.size() == freeParValues.size());
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
   assert(newFreeParams.size() == freeParValues.size());
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
      assert(fabs(sm_derivs[i] + num_deriv) < m_fracTol);

// Reset the parameters for next time around.
      srcModel.setParamValues(params_save);
      params = params_save;
   }

   std::cout << ".";
}

void LikelihoodTests::test_SourceDerivs() {
   SourceFactory * srcFactory = srcFactoryInstance();
   assert(srcFactory != 0);
   std::vector<std::string> srcNames;
   srcFactory->fetchSrcNames(srcNames);

// Create a 1 GeV Front-converting event from the anticenter at
// MET=1000 s.
   double ra(86.4);
   double dec(28.9);
   double energy(1e3);
   double time(1e3);
   ScData * scData = ScData::instance();
   assert(scData != 0);
   astro::SkyDir zAxis = scData->zAxis(time);
   double muZenith(1.);
   int type = 0;
   Event myEvent(ra, dec, energy, time, zAxis.ra(), zAxis.dec(), 
                 muZenith, type);

// Loop over Sources and check fluxDensity and Npred derivatives
// for each.
   for (int i = 0; i < 3; i++) {
      Source * src = srcFactory->create(srcNames[i]);
      assert(src != 0);
      if (src->getType() == "Diffuse") {
         DiffuseSource *diffuse_src = dynamic_cast<DiffuseSource *>(src);
         myEvent.computeResponse(*diffuse_src);
      }
      Source::FuncMap srcFuncs = src->getSrcFuncs();
      assert(srcFuncs.count("Spectrum"));

      std::vector<double> paramValues;
      srcFuncs["Spectrum"]->getFreeParamValues(paramValues);

      std::vector<std::string> paramNames;
      srcFuncs["Spectrum"]->getFreeParamNames(paramNames);
      assert(paramValues.size() == paramNames.size());

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
   std::cout << ".";
}

void LikelihoodTests::test_PointSource() {
   std::string eventFile = m_rootPath + "/data/single_src_events_0000.fits";
   std::string scDataFile = m_rootPath + "/data/single_src_scData_0000.fits";
//    std::string eventFile = m_rootPath 
//       + "/data/single_src_g25_events_0000.fits";
//    std::string scDataFile = m_rootPath 
//       + "/data/single_src_g25_scData_0000.fits";

   tearDown();
   setUp();

   std::vector<Event> events;
   readEventData(eventFile, scDataFile, events);

   SourceFactory * srcFactory = srcFactoryInstance("", scDataFile, "", "");

   Source * src = srcFactory->create("Crab Pulsar");

   RoiCuts * roiCuts = RoiCuts::instance();
   
   double Nobs = 0;
   for (unsigned int i = 0; i < events.size(); i++) {
      if (roiCuts->accept(events[i])) {
         Nobs++;
      }
   }
   std::cout << Nobs << "  "
             << src->Npred() << std::endl;

// Consider the observation over two-day intervals over the first
// twenty days, resetting the ROI accordingly and force the
// PointSource exposure to be recomputed for each interval.
   double chi2 = 0;
   double tstep = 2.*8.64e4;
   for (int j = 0; j < 10; j++) {
      double tmin = j*tstep;
      double tmax = tmin + tstep;
      RoiCuts::setCuts(86.4, 28.9, 25., 30., 3.16e5, tmin, tmax, -1.);

      src->setDir(83.57, 22.01, true, false);

      Nobs = 0;
      for (unsigned int i = 0; i < events.size(); i++) {
         if (roiCuts->accept(events[i])) {
            Nobs++;
         }
      }
      double Npred = src->Npred();
      chi2 += pow((Nobs - Npred), 2)/Nobs;
      std::cout << j << "  " 
                << Nobs << "  "
                << Npred << std::endl;
   }
   std::cout << "chi^2 = " << chi2 << std::endl;
}

void LikelihoodTests::test_DiffuseSource() {
   std::string eventFile = m_rootPath + "/data/galdiffuse_events_0000.fits";
   std::string scDataFile = m_rootPath + "/data/single_src_scData_0000.fits";
//    std::string eventFile 
//       = m_rootPath + "/data/galdiffuse_g25_events_0000.fits";
//    std::string scDataFile 
//       = m_rootPath + "/data/galdiffuse_g25_scData_0000.fits";

   std::vector<Event> events;
   readEventData(eventFile, scDataFile, events);

   astro::SkyDir anticenter(180., 0., astro::SkyDir::GALACTIC);

   double chi2 = 0;
   for (int i = 0; i < 5; i++) {
      tearDown();
      setUp();

      double tstep = 0.1*8.64e4;
      double tmin = i*tstep;
      double tmax = tmin + tstep;
      RoiCuts::setCuts(anticenter.ra(), anticenter.dec(), 20.,
                       30., 3.1623e5, tmin, tmax);

      std::ostringstream roiFile;
      roiFile << m_rootPath << "/data/RoiCuts_" << i << ".xml";

      RoiCuts * roiCuts = RoiCuts::instance();

      roiCuts->writeXml(roiFile.str(), "Anticenter region");

      std::ostringstream expMapFile;
      expMapFile << m_rootPath << "/data/expMap_" << i << ".fits";
                 
      SourceFactory * srcFactory 
         = srcFactoryInstance(roiFile.str(), scDataFile, expMapFile.str());

      Source * src = srcFactory->create("Galactic Diffuse");
      
      double Nobs = 0;
      for (unsigned int j = 0; j < events.size(); j++) {
         if (roiCuts->accept(events[j])) {
            Nobs++;
         }
      }
      double Npred = src->Npred();
      chi2 += pow((Nobs - Npred), 2)/Nobs;
      std::cout << i << "  " 
                << Nobs << "  "
                << Npred << std::endl;

      std::remove(roiFile.str().c_str());
   }
   std::cout << "chi^2 = " << chi2 << std::endl;
}

void LikelihoodTests::readEventData(const std::string &eventFile,
                                    const std::string &scDataFile,
                                    std::vector<Event> &events) {
   events.clear();
   ScData::readData(scDataFile, 2, true);
   ScData * scData = ScData::instance();

   tuple::ITable::Factory & tableFactory = *tuple::ITable::Factory::instance();
   tuple::ITable * eventTable = tableFactory(eventFile);
   int type;

   const double & ra = eventTable->selectColumn("RA");
   const double & dec = eventTable->selectColumn("DEC");
   const double & energy = eventTable->selectColumn("ENERGY");
   const double & time = eventTable->selectColumn("TIME");
   const double & zenith_angle = eventTable->selectColumn("ZENITH_ANGLE");
   const double & conversion_layer 
      = eventTable->selectColumn("CONVERSION_LAYER");
   
   for (tuple::Iterator it = eventTable->begin(); it != eventTable->end();
        ++it) {
      astro::SkyDir zAxis = scData->zAxis(time);
      if (conversion_layer < 12) {
         type = 0;
      } else {
         type = 1;
      }
      events.push_back(Event(ra, dec, energy, time, zAxis.ra(), zAxis.dec(),
                             cos(zenith_angle*M_PI/180.), type));
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
srcFactoryInstance(const std::string & roiFile,
                   const std::string & scFile,
                   const std::string & expMapFile,
                   const std::string & sourceXmlFile,
                   bool verbose) {
   if (m_srcFactory == 0) {
      RoiCuts * roiCuts = RoiCuts::instance();
      if (roiFile == "") {
         roiCuts->setCuts(m_roiFile);
      } else {
         roiCuts->setCuts(roiFile);
      }

      if (scFile == "") {
         ScData::readData(m_scFile, 2, true);
      } else {
         ScData::readData(scFile, 2, true);
      }

      if (expMapFile == "") {
         ExposureMap::readExposureFile(m_expMapFile);
      } else {
         ExposureMap::readExposureFile(expMapFile);
      }

      optimizers::FunctionFactory * funcFactory = funcFactoryInstance();

      m_srcFactory = new SourceFactory(verbose);
      if (sourceXmlFile == "") {
         m_srcFactory->readXml(m_sourceXmlFile, *funcFactory);
      } else {
         m_srcFactory->readXml(sourceXmlFile, *funcFactory);
      }
   }
   return m_srcFactory;
}      

int main() {
   LikelihoodTests unit;

   unit.test_RoiCuts();
   unit.test_SourceFactory();
   unit.test_XmlBuilders();
   unit.test_SourceModel();
   unit.test_SourceDerivs();
   unit.test_PointSource();
   unit.test_DiffuseSource();

   std::cout << "all tests ok" << std::endl;
   return 1;
}
