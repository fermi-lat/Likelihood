/** 
 * @file makeExposureCube.cxx
 * @brief Create an Exposure hypercube. This program is based on
 * map_tools/exposure_cube.cxx
 * @author J. Chiang
 *
 *  $Header$
 */

#include <sstream>
#include "facilities/Util.h"
#include "tuple/ITable.h"
#include "hoopsUtil/ParametersBase.h"
#include "hoopsUtil/RunParams.h"
#include "map_tools/ExposureHyperCube.h"

#include "Likelihood/RoiCuts.h"
#include "Likelihood/LikeExposure.h"

class ExposureMapParameters : public hoopsUtil::ParametersBase {

public:

   ExposureMapParameters(int argc, char* argv[]) 
      : hoopsUtil::ParametersBase(argc, argv) {
      m_scFile = getValue<std::string>("Spacecraft file");
      facilities::Util::expandEnvVar(&m_scFile);
      m_outputFile = getValue<std::string>("Output file");
      facilities::Util::expandEnvVar(&m_outputFile);
      getValue<double>("cos_theta step");
      getValue<double>("pixel size");
      m_roiFile = getValue<std::string>("ROI_file");
      facilities::Util::expandEnvVar(&m_roiFile);
   }

   virtual ~ExposureMapParameters() {}

   const std::string & scFile() const {return m_scFile;}
   const std::string & outputFile() const {return m_outputFile;}
   const std::string & roiFile() const {return m_roiFile;}

private:

   std::string m_scFile;
   std::string m_outputFile;
   std::string m_roiFile;
   
};

int main(int argc, char * argv[]) {
//   ExposureMapParameters pars(argc, argv);
   hoopsUtil::RunParams pars(argc, argv);

// The LikeExposure object is a data cube of exposure times in 
// (ra, dec, cost_theta):
//    Likelihood::LikeExposure exposure(pars["pixel size"], 
//                                      pars["cos_theta step"],
//                                      pars.roiFile());
   Likelihood::LikeExposure exposure(pars.getParam<double>("pixel size"), 
                                     pars.getParam<double>("cos_theta step"),
                                     pars.getParam<std::string>("ROI_file"));

// Access the Spacecraft data using as an ITable object.
   tuple::ITable::Factory & factory = *tuple::ITable::Factory::instance();
//   tuple::ITable & scData = *factory(pars.scFile(), "Ext1");
   tuple::ITable & scData 
      = *factory(pars.getParam<std::string>("Spacecraft file"), "Ext1");

// This performs the sums over the scData rows.
   exposure.load(scData);

// Create the fits output file from the Exposure file
   map_tools::ExposureHyperCube cube(exposure, 
                                     pars.getParam<std::string>("Output file"));

// Add the ROI info as HISTORY keyword data.
   Likelihood::RoiCuts * roiCuts = Likelihood::RoiCuts::instance();
   std::ostringstream roi_xml;
   roiCuts->writeXml(roi_xml);
   cube.setKey("HISTORY", roi_xml.str());

   return 0;
}
