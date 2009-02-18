/** 
 * @file SpatialMap.cxx
 * @brief Implementation of Function object class that returns interpolated
 * image values of a FITS image file.
 * 
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SpatialMap.cxx,v 1.26 2009/02/18 02:01:38 jchiang Exp $
 *
 */

#include <stdexcept>

#include "facilities/Util.h"

#include "st_stream/StreamFormatter.h"

#include "st_facilities/Util.h"

#include "Likelihood/Event.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/SkyDirArg.h"
#include "Likelihood/SpatialMap.h"
#include "Likelihood/WcsMap.h"

namespace Likelihood {

SpatialMap::SpatialMap() : m_wcsmap(0) {
   init();      
}

SpatialMap::SpatialMap(const std::string & fitsFile,
                       const std::string & extension) : m_wcsmap(0) {
   init();
   readFitsFile(fitsFile, extension);
}

SpatialMap::SpatialMap(const SpatialMap & rhs) 
   : optimizers::Function(rhs) {
   m_fitsFile = rhs.m_fitsFile;
   m_extension = rhs.m_extension;
   if (rhs.m_wcsmap) {
      m_wcsmap = new WcsMap(*(rhs.m_wcsmap));
   }
}

SpatialMap & SpatialMap::operator=(const SpatialMap & rhs) {
   if (this != &rhs) {
      optimizers::Function::operator=(rhs);
      m_fitsFile = rhs.m_fitsFile;
      m_extension = rhs.m_extension;
      delete m_wcsmap;
      if (rhs.m_wcsmap) {
         m_wcsmap = new WcsMap(*(rhs.m_wcsmap));
      }
   }
   return *this;
}

SpatialMap::~SpatialMap() {
   delete m_wcsmap;
}

void SpatialMap::init() {
// This Function has one Parameter, an overall normalization, 
// but set it to be unit constant.
   int nParams = 1;
   setMaxNumParams(nParams);
   m_genericName = "SpatialMap";
   addParam("Prefactor", 1, false);
   setParamAlwaysFixed("Prefactor");
   m_normParName = "Prefactor";
}

void SpatialMap::readFitsFile(const std::string & fitsFile,
                              const std::string & extension) {
   m_fitsFile = fitsFile;
   m_extension = extension;

   std::string expandedFileName(fitsFile);

   facilities::Util::expandEnvVar(&expandedFileName);

   if (!st_facilities::Util::fileExists(expandedFileName)) {
// The following to StreamFormatter is necessary since Xerces seems to
// corrupt the exception handling when this method is called from
// SourceFactory::readXml and the program simply aborts.
      st_stream::StreamFormatter formatter("SpatialMap", "readFitsFile", 2);
      formatter.err() << "File not found: " << expandedFileName << std::endl;
      throw std::runtime_error("File not found: " + expandedFileName);
   }
   m_wcsmap = new WcsMap(expandedFileName, extension);
}

double SpatialMap::value(optimizers::Arg & arg) const {
   astro::SkyDir dir;
   dynamic_cast<SkyDirArg &>(arg).fetchValue(dir);
   return value(dir);
}

double SpatialMap::value(const astro::SkyDir & dir) const {
   return m_wcsmap->operator()(dir);
}

double SpatialMap::diffuseResponse(const ResponseFunctions & respFuncs,
                                   const Event & event) const {
   double my_value(0);
   for (int ilon(1); ilon < m_wcsmap->nxpix()+1; ilon++) {
      for (int ilat(1); ilat < m_wcsmap->nypix()+1; ilat++) {
         astro::SkyDir srcDir(m_wcsmap->skyDir(ilon, ilat));
         double totalResponse = 
            respFuncs.totalResponse(event.getEnergy(), event.getEnergy(),
                                    event.zAxis(), event.xAxis(),
                                    srcDir, event.getDir(), event.getType());
         double addend(m_wcsmap->pixelValue(ilon, ilat)*
                       m_wcsmap->solidAngle(ilon, ilat)*
                       totalResponse);
         my_value += addend;
//          std::cout << ilon << "  "
//                    << ilat << "  "
//                    << srcDir.ra() << "  "
//                    << srcDir.dec() << "  "
//                    << m_wcsmap->pixelValue(ilon, ilat) << "  "
//                    << m_wcsmap->solidAngle(ilon, ilat) << "  "
//                    << totalResponse << "  "
//                    << addend << std::endl;
      }
   }
   return my_value;
}

bool SpatialMap::insideMap(const astro::SkyDir & dir) const {
   return m_wcsmap->insideMap(dir);
}

std::pair<astro::SkyDir, astro::SkyDir> 
SpatialMap::minMaxDistPixels(const astro::SkyDir & dir) const {
   return m_wcsmap->minMaxDistPixels(dir);
}

void SpatialMap::getCorners(std::vector<astro::SkyDir> & corners) const {
   m_wcsmap->getCorners(corners);
}

} // namespace Likelihood
