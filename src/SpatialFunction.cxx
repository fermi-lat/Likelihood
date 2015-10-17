/** 
 * @file SpatialFunction.cxx
 * @brief Implementation of Function object class that represents an
 * analytic function in 2 spatial dimensions.  This class serves as
 * a base class for various spatial template classes.
 * 
 * @author M. Wood
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SpatialFunction.cxx,v 1.37 2015/03/21 05:38:04 jchiang Exp $
 *
 */

#include <algorithm>
#include <iostream>
#include <stdexcept>

#include "facilities/Util.h"

#include "st_stream/StreamFormatter.h"

#include "st_facilities/Util.h"

#include "astro/SkyDir.h"

#include "optimizers/dArg.h"
#include "optimizers/Function.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/Event.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/SpatialFunction.h"
#include "Likelihood/SkyDirArg.h"

namespace Likelihood {

SpatialFunction::SpatialFunction(const std::string& name,unsigned nparam) 
  : optimizers::Function(name, nparam, "") {
  m_dir = astro::SkyDir(0.0, 0.0);
  init(0.0,0.0);
}

SpatialFunction::SpatialFunction(const std::string& name,unsigned nparam, 
				 double ra, double dec) 
  : optimizers::Function(name, nparam, "") {
  m_dir = astro::SkyDir(ra, dec);
  init(ra,dec);
}

void SpatialFunction::init(double ra, double dec) {
  addParam("RA", ra, false);
  parameter("RA").setBounds(-360., 360.);
  setParamAlwaysFixed("RA");
  addParam("DEC", dec, false);
  parameter("DEC").setBounds(-90., 90.);
  setParamAlwaysFixed("DEC");
}

void SpatialFunction::setCenter(double ra, double dec) {
  m_dir = astro::SkyDir(ra, dec);
}

SpatialFunction::SpatialFunction(const SpatialFunction & rhs) 
  : optimizers::Function(rhs), m_dir(rhs.m_dir) {
}

SpatialFunction & 
SpatialFunction::operator=(const SpatialFunction & rhs) {
   if (this != &rhs) {
      optimizers::Function::operator=(rhs);
      m_dir = rhs.m_dir;
   }
   return *this;
}

SpatialFunction::~SpatialFunction() {
}

const astro::SkyDir& SpatialFunction::dir() const {
  return m_dir;
}

void SpatialFunction::setParams(const DOMElement * elt) {
  optimizers::Function::setParams(elt);
  update();
}

void SpatialFunction::update() {
  double ra = getParam("RA").getValue();
  double dec = getParam("DEC").getValue();  
  setCenter(ra,dec);
}

} // namespace Likelihood
