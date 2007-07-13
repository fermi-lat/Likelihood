/** 
 * @file SpatialMap.cxx
 * @brief Implementation of Function object class that returns interpolated
 * image values of a FITS image file.
 * 
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SpatialMap.cxx,v 1.21 2007/03/14 20:11:10 jchiang Exp $
 *
 */

#include <stdexcept>

#include "facilities/Util.h"

#include "st_stream/StreamFormatter.h"

#include "st_facilities/Util.h"

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

} // namespace Likelihood
