/** 
 * @file SpatialMap.cxx
 * @brief Implementation of Function object class that returns interpolated
 * image values of a FITS image file.
 * 
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SpatialMap.cxx,v 1.17 2005/10/04 05:38:08 jchiang Exp $
 *
 */

#include "facilities/Util.h"

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
}

void SpatialMap::readFitsFile(const std::string & fitsFile,
                              const std::string & extension) {
   m_fitsFile = fitsFile;
   m_extension = extension;

   facilities::Util::expandEnvVar(&m_fitsFile);

   m_wcsmap = new WcsMap(m_fitsFile, extension);
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
