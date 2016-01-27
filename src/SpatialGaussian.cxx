/** 
 * @file SpatialGaussian.cxx
 * @brief Implementation of Function object class for a 2D spatial gaussian.
 * 
 * @author M. Wood
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/SpatialGaussian.cxx,v 1.2 2015/12/03 23:47:46 mdwood Exp $
 *
 */

#include <algorithm>
#include <iostream>
#include <stdexcept>

#include "gsl/gsl_sf_bessel.h"

#include "facilities/Util.h"

#include "st_stream/StreamFormatter.h"

#include "st_facilities/GaussianQuadrature.h"
#include "st_facilities/Util.h"

#include "astro/SkyDir.h"

#include "optimizers/dArg.h"
#include "optimizers/Function.h"

#include "Likelihood/Event.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/SpatialGaussian.h"
#include "Likelihood/SpatialFunction.h"
#include "Likelihood/SkyDirArg.h"

namespace {

  double gauss(double x,double sigma) {

    return 1./(2*M_PI*sigma*sigma)*std::exp(-x*x/(2*sigma*sigma));

   }
}

namespace Likelihood {

int SpatialGaussian::RadialIntegrand::s_ncall = 0;

double SpatialGaussian::RadialIntegrand::operator()(double xp) const {

  const double s2 = std::pow(m_sigma,2);
  double xx = m_x*xp/s2;
  double je = gsl_sf_bessel_I0_scaled(xx);      

  s_ncall += 1;

  return xp*m_fn(m_energy,xp)*je*std::exp(xx - (m_x*m_x+xp*xp)/(2*s2))/s2;
}

double SpatialGaussian::convolve(const ResponseFunctor& fn, 
				 double energy, double x, 
				 double sigma, double err) {

  const double xmin = std::max(x - 6*sigma,0.0);
  const double xmax = x + 6*sigma;
  RadialIntegrand rIntegrand(fn,energy,x,sigma);

  int ierr;
  double v0 = 0, v1 = 0;

  if(xmin < 0)
    v0 = st_facilities::GaussianQuadrature::dgaus8(rIntegrand, xmin,
						   0, err, ierr);
  v1 = st_facilities::GaussianQuadrature::dgaus8(rIntegrand, 0,
						 xmax, err, ierr);
  
  return v0 + v1;
}

SpatialGaussian::SpatialGaussian() 
  : SpatialFunction("SpatialGaussian",3) {
  m_sigma = 1.0;
  addParam("Sigma", 1.0, false);
  parameter("Sigma").setBounds(0.0, 180.);
}

SpatialGaussian::SpatialGaussian(double ra, double dec, double sigma) 
  : SpatialFunction("SpatialGaussian",3,ra,dec) {
  m_sigma = sigma;
  addParam("Sigma", m_sigma, false);
  parameter("Sigma").setBounds(0.0, 180.);
}

SpatialGaussian::SpatialGaussian(const SpatialGaussian & rhs) 
  : SpatialFunction(rhs), m_sigma(rhs.m_sigma) {

}

SpatialGaussian & SpatialGaussian::operator=(const SpatialGaussian & rhs) {
   if (this != &rhs) {
      SpatialFunction::operator=(rhs);
      m_sigma = rhs.m_sigma;
   }
   return *this;
}

SpatialGaussian::~SpatialGaussian() {
}

double SpatialGaussian::value(const astro::SkyDir & dir) const {
   double separation = this->dir().difference(dir)*180./M_PI;
   return gauss(separation,m_sigma)*std::pow(M_PI/180.,-2);
}

double SpatialGaussian::value(double separation, double sigma) const {
   return gauss(separation,sigma)*std::pow(M_PI/180.,-2);
}

double SpatialGaussian::spatialResponse(const astro::SkyDir & dir, double energy, const MeanPsf& psf) 
  const {
   double separation = dir.difference(this->dir())*180./M_PI;
   return SpatialGaussian::convolve(BinnedResponseFunctor(psf),energy,separation,m_sigma);
}

double SpatialGaussian::spatialResponse(double separation, double energy, const MeanPsf& psf) const {
  return SpatialGaussian::convolve(BinnedResponseFunctor(psf),energy,separation,m_sigma);
}

double SpatialGaussian::diffuseResponse(const ResponseFunctor& fn, double energy,
					double separation) const {
  return convolve(fn,energy,separation,m_sigma);
}

double SpatialGaussian::getDiffRespLimits(const astro::SkyDir & dir, 
					  double & mumin, double & mumax,
					  double & phimin, double & phimax) const {
   mumin = std::cos(2.*dir.difference(this->dir()) + 3.*m_sigma*M_PI/180.);
   mumax = 1;
   phimin = 0;
   phimax = 2.*M_PI;   
}

void SpatialGaussian::update() {
   SpatialFunction::update();
   double sigma = getParam("Sigma").getValue();  
   m_sigma = sigma;
}

double SpatialGaussian::value(const optimizers::Arg & x) const {

   const SkyDirArg & dir = dynamic_cast<const SkyDirArg &>(x);
   double offset = dir().difference(this->dir())*180./M_PI;
   return value(offset,m_sigma);
}

double SpatialGaussian::derivByParamImp(const optimizers::Arg & x, 
                                      const std::string & parName) const {

   std::ostringstream message;
   message << "SpatialGaussian: cannot take derivative wrt "
           << "parameter " << parName;
   throw std::runtime_error(message.str());
}

} // namespace Likelihood
