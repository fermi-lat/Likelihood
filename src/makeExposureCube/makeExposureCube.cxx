/** 
 * @file makeExposureCube.cxx
 * @brief Create an Exposure hypercube. This program is based on
 * map_tools/exposure_cube.cxx
 * @author J. Chiang
 *
 *  $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/test/makeExposureCube.cxx,v 1.3 2004/03/19 04:56:15 jchiang Exp $
 */

#include <sstream>

#include "facilities/Util.h"

#include "tuple/ITable.h"

#include "hoops/hoops.h"
#include "hoops/hoops_prompt_group.h"

#include "map_tools/ExposureHyperCube.h"

#include "Likelihood/RoiCuts.h"
#include "Likelihood/LikeExposure.h"

int main(int argc, char * argv[]) {
   hoops::ParPromptGroup pars(argc, argv);
   pars.Prompt();
   pars.Save();

// The LikeExposure object is a data cube of exposure times in 
// (ra, dec, cost_theta):
   Likelihood::LikeExposure exposure(pars["pixel size"], 
                                     pars["cos_theta step"], 
                                     pars["ROI_file"]);

// Access the Spacecraft data using as an ITable object.
   tuple::ITable::Factory & factory = *tuple::ITable::Factory::instance();
   tuple::ITable & scData = *factory(pars["Spacecraft file"], "Ext1");

// This performs the sums over the scData rows.
   exposure.load(scData);

// Create the fits output file from the Exposure file
   map_tools::ExposureHyperCube cube(exposure, pars["Output file"]);
   
// Add the ROI info as HISTORY keyword data.
   Likelihood::RoiCuts * roiCuts = Likelihood::RoiCuts::instance();
   std::ostringstream roi_xml;
   roiCuts->writeXml(roi_xml);
   cube.setKey("HISTORY", roi_xml.str());

   return 0;
}
