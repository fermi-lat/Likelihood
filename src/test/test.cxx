/**
 * @file test.cxx
 * @brief Test program for Likelihood.  Use CppUnit-like idioms.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/test/test.cxx,v 1.1 2004/02/20 00:02:08 jchiang Exp $
 */

#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <cmath>
#include <cassert>

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

#include "optimizers/FunctionFactory.h"

#include "latResponse/AcceptanceCone.h"
#include "latResponse/IrfsFactory.h"

#include "Likelihood/ExposureMap.h"
#include "Likelihood/FluxBuilder.h"
#include "Likelihood/SourceModelBuilder.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/Source.h"
#include "Likelihood/SourceFactory.h"
#include "Likelihood/SourceModel.h"
#include "Likelihood/SpatialMap.h"

using namespace Likelihood;
using optimizers::Parameter;

class SourceData {

public:

   SourceData() {
      char * srcNames[] = {"Extragalactic Diffuse", "Galactic Diffuse",
                           "PKS 0528+134", "Crab Pulsar", "Geminga"};
      char * srcTypes[] = {"Diffuse", "Diffuse", "Point", "Point", "Point"};
      char * spatialModels[] = {"ConstantValue", "SpatialMap", 
                                "SkyDirFunction", "SkyDirFunction", 
                                "SkyDirFunction"};
      for (unsigned int i = 0; i < 5; i++) {
         m_srcId[srcNames[i]] = i;
         m_srcTypes.push_back(srcTypes[i]);
         m_spatialModels.push_back(spatialModels[i]);
      }
      setParameters();
   }

   ~SourceData() {}

   void getSrcNames(std::vector<std::string> &srcNames) {
      srcNames.clear();
      for (std::map<std::string, unsigned int>::iterator it = m_srcId.begin();
           it != m_srcId.end(); it++) {
         srcNames.push_back(it->first);
      }
   }

   const std::string & srcType(const std::string & srcName) {
      assert(m_srcId.count(srcName));
      return m_srcTypes[m_srcId[srcName]];
   }

   const std::string & spatialModel(const std::string & srcName) {
      assert(m_srcId.count(srcName));
      return m_spatialModels[m_srcId[srcName]];
   }

   const std::vector<Parameter> & parameters(const std::string & srcName) {
      assert(m_srcId.count(srcName));
      return m_parameters[m_srcId[srcName]];
   }

   const Parameter & paramObject(const std::string & srcName,
                                 const std::string & paramName) {
      assert(m_srcId.count(srcName));
      const std::vector<Parameter> & params = parameters(srcName);
      for (unsigned int i = 0; i < params.size(); i++) {
         if (params[i].getName() == paramName) {
            return params[i];
         }
      }
      bool paramFound(false);
      assert(paramFound);
      return params[0];
   }

private:

   std::map<std::string, unsigned int> m_srcId;
   std::vector<std::string> m_srcTypes;
   std::vector<std::string> m_spatialModels;

   std::vector< std::vector<Parameter> > m_parameters;
   void setParameters();

};

void SourceData::setParameters() {
// We need to do this by hand to ensure an unbiased test of the source
// data. Here we accumulate parameter data only for the spectral model
// components.

   m_parameters.resize(5);
   unsigned int id = m_srcId["Extragalactic Diffuse"];

   m_parameters[id].push_back(Parameter("Prefactor", 1.32, 1e-5, 1e2, true));
   m_parameters[id].push_back(Parameter("Index", -2.1, -1., -3.5, false));
   m_parameters[id].push_back(Parameter("Scale", 1e2, 50., 2e2, false));

   id = m_srcId["Galactic Diffuse"];
   m_parameters[id].push_back(Parameter("Prefactor", 11., 1e-3, 1e3, true));
   m_parameters[id].push_back(Parameter("Index", -2.1, -1., -3.5, false));
   m_parameters[id].push_back(Parameter("Scale", 1e2, 50., 2e2, false));

   id = m_srcId["PKS 0528+134"];
   m_parameters[id].push_back(Parameter("Prefactor", 13.65, 1e-3, 1e3, true));
   m_parameters[id].push_back(Parameter("Index", -2.46, -1., -3.5, true));
   m_parameters[id].push_back(Parameter("Scale", 1e2, 30., 2e3, false));

   id = m_srcId["Crab Pulsar"];
   m_parameters[id].push_back(Parameter("Prefactor", 27., 1e-3, 1e3, true));
   m_parameters[id].push_back(Parameter("Index", -2.19, -1., -3.5, true));
   m_parameters[id].push_back(Parameter("Scale", 1e2, 30., 2e3, false));

   id = m_srcId["Geminga"];
   m_parameters[id].push_back(Parameter("Prefactor", 23.29, 1e-3, 1e3, true));
   m_parameters[id].push_back(Parameter("Index", -1.66, -1., -3.5, true));
   m_parameters[id].push_back(Parameter("Scale", 1e2, 30., 2e3, false));
}

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

private:

   std::string m_rootPath;
   double m_fracTol;

// File names for m_srcFactory.
   std::string m_roiFile;
   std::string m_expMapFile;
   std::string m_sourceXmlFile;

// Test data for m_srcFactory.
   SourceData m_srcData;

   optimizers::FunctionFactory * m_funcFactory;
   optimizers::FunctionFactory * funcFactoryInstance();

   SourceFactory * m_srcFactory;
   SourceFactory * srcFactoryInstance();
   
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
   
// Fractional tolerance for double comparisons.
   m_fracTol = 1e-4;

// Use lazy evaluation for m_funcFactory and m_srcFactory.
   m_funcFactory = 0;
   m_srcFactory = 0;

   m_roiFile = m_rootPath + "/data/anticenter_Roi.xml";
   m_expMapFile = m_rootPath + "/data/anticenter_expMap.fits";
   m_sourceXmlFile = m_rootPath + "/data/anticenter_model.xml";
}

void LikelihoodTests::tearDown() {
   delete m_funcFactory;
   delete m_srcFactory;
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
      ASSERT_EQUALS(srcFuncs["Spectrum"]->getParamValue("Prefactor"),
                    m_srcData.paramObject(name, "Prefactor").getValue());
      ASSERT_EQUALS(srcFuncs["Spectrum"]->getParamValue("Index"), 
                    m_srcData.paramObject(name, "Index").getValue());
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
   std::string fluxXmlFile = m_rootPath + "/data/fluxBuilder.xml";
   fluxBuilder.write(fluxXmlFile);

   std::string srcModelXmlFile = m_rootPath + "/data/srcModelBuilder.xml";
   srcModelBuilder.write(srcModelXmlFile);

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
   for (unsigned int i = 0; i < srcNames.size(); i++) {
      std::string name = srcNames[i];
      ASSERT_EQUALS(srcModel.getParamValue("Prefactor", "Spectrum", name),
                    m_srcData.paramObject(name, "Prefactor").getValue());
      ASSERT_EQUALS(srcModel.getParamValue("Index", "Spectrum", name),
                    m_srcData.paramObject(name, "Index").getValue());
      ASSERT_EQUALS(srcModel.getParamValue("Scale", "Spectrum", name),
                    m_srcData.paramObject(name, "Scale").getValue());
   }

   std::cout << ".";
}

optimizers::FunctionFactory * LikelihoodTests::funcFactoryInstance() {
   if (m_funcFactory == 0) {
      m_funcFactory = new optimizers::FunctionFactory();
      m_funcFactory->addFunc("SkyDirFunction", new SkyDirFunction(), false);
      m_funcFactory->addFunc("SpatialMap", new SpatialMap(), false);
   }
   return m_funcFactory;
}

SourceFactory * LikelihoodTests::srcFactoryInstance() {
   if (m_srcFactory == 0) {
      RoiCuts * roiCuts = RoiCuts::instance();
      roiCuts->setCuts(m_roiFile);

      ExposureMap::readExposureFile(m_expMapFile);

      optimizers::FunctionFactory * funcFactory = funcFactoryInstance();

      m_srcFactory = new SourceFactory();
      m_srcFactory->readXml(m_sourceXmlFile, *funcFactory);
   }
   return m_srcFactory;
}      

int main() {
   LikelihoodTests unit;

   unit.test_RoiCuts();
   unit.test_SourceFactory();
   unit.test_XmlBuilders();
   unit.test_SourceModel();

   std::cout << "all tests ok" << std::endl;
   return 1;
}
