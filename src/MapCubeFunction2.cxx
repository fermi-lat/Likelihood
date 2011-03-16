/**
 * @file MapCubeFunction2.cxx
 * @brief Encapsulation of 3D FITS image of a diffuse source with 
 * position-dependent spectral variation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/MapCubeFunction2.cxx,v 1.34 2011/02/09 23:11:50 jchiang Exp $
 */

#include <cmath>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "facilities/Util.h"

#include "st_stream/StreamFormatter.h"

#include "Likelihood/MapCubeFunction2.h"
#include "Likelihood/SkyDirArg.h"
#include "Likelihood/WcsMap2.h"

namespace Likelihood {

MapCubeFunction2::MapCubeFunction2() 
   : m_wcsmap(0) {
   init();
}

MapCubeFunction2::MapCubeFunction2(const std::string & fitsFile) 
   : optimizers::Function(), MapBase(fitsFile), 
     m_wcsmap(new WcsMap2(fitsFile)) {
   init();
}

MapCubeFunction2::MapCubeFunction2(const MapCubeFunction2 & rhs)
   : optimizers::Function(rhs), MapBase(rhs), 
     m_wcsmap(new WcsMap2(*rhs.m_wcsmap)) {
   init();
}

MapCubeFunction2 & MapCubeFunction2::operator=(const MapCubeFunction2 & rhs) {
   if (this != &rhs) {
      optimizers::Function::operator=(rhs);
      MapBase::operator=(rhs);
      m_wcsmap = new WcsMap2(*rhs.m_wcsmap);
   }
   return *this;
}

MapCubeFunction2::~MapCubeFunction2() {
   delete m_wcsmap;
}

double MapCubeFunction2::value(optimizers::Arg & xarg) const {
   SkyDirArg & dir = dynamic_cast<SkyDirArg &>(xarg);
   double energy = dir.energy();
   double value = m_wcsmap->operator()(dir(), energy);
   return value*getParam("Normalization").getTrueValue();
}

void MapCubeFunction2::init() {
   setMaxNumParams(1);
// Leave this parameter fixed, modifying the overall normalization
// via a ConstantValue function as the spectral component.
   addParam("Normalization", 1, false);

   m_funcType = Addend;
   m_argType = "";
   m_genericName = "MapCubeFunction2";
   m_normParName = "Normalization";
}


double MapCubeFunction2::mapIntegral() const {
   return m_wcsmap->mapIntegral();
}

}
