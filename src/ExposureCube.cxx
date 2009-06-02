/**
 * @file ExposureCube.cxx
 * @brief Implementation for ExposureCube wrapper class of map_tools::Exposure
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/ExposureCube.cxx,v 1.5 2009/06/02 06:21:43 jchiang Exp $
 */

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "Likelihood/ExposureCube.h"

namespace Likelihood {

ExposureCube::ExposureCube(const ExposureCube & other) 
   : m_exposure(new map_tools::Exposure(*(other.m_exposure))),
     m_weightedExposure(0), m_haveFile(other.m_haveFile),
     m_fileName(other.m_fileName), 
     m_hasPhiDependence(other.m_hasPhiDependence) {
   if (other.m_weightedExposure) {
      m_weightedExposure = new map_tools::Exposure(*(other.m_weightedExposure));
   }
}

void ExposureCube::readExposureCube(std::string filename) {
   facilities::Util::expandEnvVar(&filename);
   m_fileName = filename;
   m_exposure = new map_tools::Exposure(filename);
   try {
      m_weightedExposure = new map_tools::Exposure(filename,
                                                   "WEIGHTED_EXPOSURE");
   } catch(tip::TipException &) {
      m_weightedExposure = 0;
   }
   m_haveFile = true;
   m_hasPhiDependence = phiDependence(filename);
}

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
