/**
 * @file test.cxx
 * @brief Test program for Likelihood.  Use CppUnit-like idioms.
 * @author J. Chiang
 * 
 * $Header$
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
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/SourceModel.h"
#include "Likelihood/SpatialMap.h"

using namespace Likelihood;

class LikelihoodTests {

public:

   LikelihoodTests() {
      setUp();
   }

   ~LikelihoodTests() {
      tearDown();
   }

   void setUp() {
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
   }

   void tearDown() {
   }

   void test_RoiCuts();
   void test_XmlBuilders();

private:

   std::string m_rootPath;
   double m_fracTol;
   
};

#define ASSERT_EQUALS(X, Y) assert(fabs( (X - Y)/Y ) < m_fracTol)

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

void LikelihoodTests::test_XmlBuilders() {

   SourceModel srcModel;
   
// Create the FunctionFactory and SourceFactory.
   optimizers::FunctionFactory funcFactory;
   try {
// Add the Functions needed for spatial modeling.
      funcFactory.addFunc("SkyDirFunction", new SkyDirFunction(), false);
      funcFactory.addFunc("SpatialMap", new SpatialMap(), false);
   } catch (optimizers::Exception &eObj) {
      std::cout << eObj.what() << std::endl;
   }

   std::string expFile = m_rootPath + "/data/anticenter_expMap.fits";
   ExposureMap::readExposureFile(expFile);

   std::string modelFile = m_rootPath + "/data/anticenter_model.xml";
   srcModel.readXml(modelFile, funcFactory);

   std::string newModelFile = m_rootPath + "/data/newSourceModel.xml";
   srcModel.writeXml(newModelFile);
   srcModel.reReadXml(newModelFile);
   srcModel.deleteAllSources();
   srcModel.readXml(newModelFile, funcFactory);

   std::string fluxFile = m_rootPath + "/data/flux_model.xml";
   srcModel.write_fluxXml(fluxFile);

   std::cout << ".";
}

int main() {
   LikelihoodTests unit;

   unit.test_RoiCuts();
   unit.test_XmlBuilders();

   std::cout << std::endl;
   return 1;
}
