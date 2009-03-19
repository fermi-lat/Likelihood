/**
 * @file ExposureCube.cxx
 * @brief Implementation for ExposureCube wrapper class of map_tools::Exposure
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/ExposureCube.cxx,v 1.3 2009/03/16 23:19:07 jchiang Exp $
 */

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "Likelihood/ExposureCube.h"

namespace Likelihood {

bool ExposureCube::phiDependence(const std::string & filename) const {
   const tip::Table * table 
      = tip::IFileSvc::instance().readTable(filename, "EXPOSURE");
   const tip::Header & header(table->getHeader());
   long nphibins;
   try {
      header["PHIBINS"].get(nphibins);
   } catch (tip::TipException &) {
      nphibins = 0;
   }
   return nphibins > 0;
}

} // namespace Likelihood
