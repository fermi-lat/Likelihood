/**
 * @file test.cxx
 * @brief Test program for Likelihood.  Use CppUnit-like idioms.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/test/test.cxx,v 1.4 2004/02/21 05:11:52 jchiang Exp $
 */

#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <cmath>
#include <cassert>
#include <cstdio>

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

#include "facilities/Util.h"

#include "optimizers/FunctionFactory.h"

#include "latResponse/AcceptanceCone.h"
#include "latResponse/IrfsFactory.h"

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

using namespace Likelihood;
using optimizers::Parameter;

namespace {
   void makeDomElementMap(const DomElement & rootElt,
                          std::map<std::string, DomElement> & domMap) {
      std::vector<DomElement> elts;
      xml::Dom::getChildrenByTagName(rootElt, "source", elts);
      domMap.clear();
      for (unsigned int i = 0; i < elts.size(); i++) {
         std::string name = xml::Dom::getAttribute(elts[i], "name");
         domMap[name] = elts[i];
      }
   }
   void readlines(std::string inputFile, std::vector<std::string> lines) {
      facilities::Util::expandEnvVar(&inputFile);
      std::ifstream file(inputFile.c_str());
      lines.clear();
      std::string line;
      while (std::getline(file, line, '\n')) {
         lines.push_back(line);
      }
   }
   bool compareXmlFiles(const std::string &file1,
                        const std::string &file2) {

      xml::XmlParser * parser = new xml::XmlParser();
      DomDocument doc1 = parser->parse(file1.c_str());
      DomDocument doc2 = parser->parse(file2.c_str());

// Re-serialize both docs as temporary files for line-by-line
// comparison.  We need to serialize the children of the root element
// explicitly to ensure they are written in the same order.
      std::map<std::string, DomElement> domMap1;
      makeDomElementMap(doc1.getDocumentElement(), domMap1);
      std::map<std::string, DomElement> domMap2;
      makeDomElementMap(doc2.getDocumentElement(), domMap2);

      std::string file1name = file1 + "_reserialized";
      std::ofstream firstFile(file1name.c_str());
      std::string file2name = file2 + "_reserialized";
      std::ofstream secondFile(file2name.c_str());

      std::map<std::string, DomElement>::iterator it;
      for (it = domMap1.begin(); it != domMap1.end(); it++) {
         std::string name = it->first;
         DomNode node1 = it->second.cloneNode(true);
         xml::Dom::prettyPrintElement(node1, firstFile, std::string(""));
         assert(domMap2.count(name));
         DomNode node2 = domMap2[name].cloneNode(true);
         xml::Dom::prettyPrintElement(node2, secondFile, std::string(""));
      }
      firstFile.close();
      secondFile.close();
      std::vector<std::string> file1_lines;
      readlines(file1name, file1_lines);
      std::vector<std::string> file2_lines;
      readlines(file2name, file2_lines);
      assert(file1_lines.size() == file2_lines.size());
      for (unsigned int i = 0; i < file1_lines.size(); i++) {
         if (file1_lines[i] != file2_lines[i]) return false;
      }
      std::remove(file1name.c_str());
      std::remove(file2name.c_str());
      return true;
   }
}

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

   const std::vector<std::string> & paramNames() const {
      return m_paramNames;
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

   std::vector<std::string> m_paramNames;
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

   m_paramNames.push_back("Prefactor");
   m_paramNames.push_back("Index");
   m_paramNames.push_back("Scale");
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
   void test_PointSource();

private:

   std::string m_rootPath;
   double m_fracTol;

// File names for m_srcFactory.
   std::string m_roiFile;
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

   assert(::compareXmlFiles(m_sourceXmlFile, m_srcModelXmlFile));

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

   std::cout << ".";
}

void LikelihoodTests::test_PointSource() {
   SourceFactory * srcFactory = srcFactoryInstance();
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

      std::string scFile = m_rootPath + "/data/oneday_scData_0000.fits";
      ScData::readData(scFile, 2, true);

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
   unit.test_PointSource();

   std::cout << "all tests ok" << std::endl;
   return 1;
}
