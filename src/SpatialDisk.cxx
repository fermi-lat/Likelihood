/** 
 * @file SpatialDisk.cxx
 * @brief Implementation of Function object class to represent a 2D spatial disk.
 * 
 * @author M. Wood
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SpatialDisk.cxx,v 1.37 2015/03/21 05:38:04 jchiang Exp $
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

#include "Likelihood/Event.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/SpatialDisk.h"
#include "Likelihood/SkyDirArg.h"

namespace {

  double disk(double x,double sigma) {
    if(x/sigma < 1.0) 
      return std::pow(sigma*sigma*M_PI,-1);
    else 
      return 0.0;
  }
}

namespace Likelihood {

double SpatialDisk::convolve(const Likelihood::MeanPsf& psf, 
			     double energy, double x, 
			     double sigma, int n) {

  const double xpmin = std::max(x - sigma,0.0);
  const double xpmax = x + sigma;
  const double dxp = (xpmax-xpmin)/double(n+1);
  const double s2 = sigma*sigma;

  double s = 0;
  for(int i = 0; i < n; i++)
    {
      double xp = xpmin + (double(i)+0.5)*dxp;
      double dphi = 2.0*M_PI;
      if((xp+x)/sigma>1.0 && x > 0.0)
	dphi = 2.0*std::acos((std::pow(x,2)+std::pow(xp,2)-s2)/(2*x*xp));
      s += xp*psf(energy,xp)*dphi;
    }
  
  return s/(M_PI*s2)*dxp;
}

SpatialDisk::SpatialDisk() : SpatialFunction("SpatialDisk",3) {
  m_width = 1.0;
  addParam("Width", 1.0, false);
  parameter("Width").setBounds(0.0, 180.);
}

SpatialDisk::SpatialDisk(double ra, double dec, double width) 
  : SpatialFunction("SpatialDisk",3,ra,dec) {
  m_width = width;
  addParam("Width", m_width, false);
  parameter("Width").setBounds(0.0, 180.);
}

SpatialDisk::SpatialDisk(const SpatialDisk & rhs) 
  : SpatialFunction(rhs), m_width(rhs.m_width) {
}

SpatialDisk & SpatialDisk::operator=(const SpatialDisk & rhs) {
   if (this != &rhs) {
      SpatialFunction::operator=(rhs);
      m_width = rhs.m_width;
   }
   return *this;
}

SpatialDisk::~SpatialDisk() {
}

double SpatialDisk::value(const astro::SkyDir & dir) const {
   double delta = this->dir().difference(dir)*180./M_PI;
   return disk(delta,m_width)*std::pow(M_PI/180.,2);
}

double SpatialDisk::value(double delta, double width) const {
   return disk(delta,width)*std::pow(M_PI/180.,2);
}

double SpatialDisk::value(const astro::SkyDir & dir, double energy, const MeanPsf& psf) const {
  double delta = dir.difference(this->dir())*180./M_PI;
  return SpatialDisk::convolve(psf,energy,delta,m_width);
}

double SpatialDisk::value(double delta, double energy, const MeanPsf& psf) const {
  return SpatialDisk::convolve(psf,energy,delta,m_width);
}

void SpatialDisk::update() {
  SpatialFunction::update();
  double width = getParam("Width").getValue();  
  m_width = width;
}

double SpatialDisk::value(const optimizers::Arg & x) const {
   const SkyDirArg & dir = dynamic_cast<const SkyDirArg &>(x);
   double offset = dir().difference(this->dir())*180./M_PI;
   return value(offset,m_width);
}

double SpatialDisk::derivByParamImp(const optimizers::Arg & x, 
                                      const std::string & parName) const {

   std::ostringstream message;
   message << "SpatialDisk: cannot take derivative wrt "
           << "parameter " << parName;
   throw std::runtime_error(message.str());
}

} // namespace Likelihood
