/** 
 * @file SpatialMap.cxx
 * @brief Implementation of Function object class that returns interpolated
 * image values of a FITS image file.
 * 
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SpatialMap.cxx,v 1.27 2009/02/18 06:57:43 jchiang Exp $
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

SpatialMap::SpatialMap() : optimizers::Function(), MapBase() {
   init();      
}

SpatialMap::SpatialMap(const std::string & fitsFile,
                       const std::string & extension) 
   : optimizers::Function(), MapBase(fitsFile, extension) {
   init();
}

SpatialMap::SpatialMap(const SpatialMap & rhs) 
   : optimizers::Function(rhs), MapBase(rhs) {
}

SpatialMap & SpatialMap::operator=(const SpatialMap & rhs) {
   if (this != &rhs) {
      optimizers::Function::operator=(rhs);
      MapBase::operator=(rhs);
   }
   return *this;
}

SpatialMap::~SpatialMap() {
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

double SpatialMap::value(optimizers::Arg & arg) const {
   astro::SkyDir dir;
   dynamic_cast<SkyDirArg &>(arg).fetchValue(dir);
   return value(dir);
}

double SpatialMap::value(const astro::SkyDir & dir) const {
   return m_wcsmap->operator()(dir);
}

} // namespace Likelihood
