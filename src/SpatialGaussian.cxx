/** 
 * @file SpatialGaussian.cxx
 * @brief Implementation of Function object class for a 2D spatial gaussian.
 * 
 * @author M. Wood
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SpatialGaussian.cxx,v 1.37 2015/03/21 05:38:04 jchiang Exp $
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
#include "Likelihood/SpatialGaussian.h"
#include "Likelihood/SkyDirArg.h"

#define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-08

namespace {

  struct cheb_series {
    double * c;   /* coefficients                */
    int order;    /* order of expansion          */
    double a;     /* lower interval point        */
    double b;     /* upper interval point        */
    int order_sp; /* effective single precision order */
  };
  //  typedef struct cheb_series_struct cheb_series;

  static double bi0_data[12] = {
    -.07660547252839144951,
    1.92733795399380827000,
    .22826445869203013390, 
    .01304891466707290428,
    .00043442709008164874,
    .00000942265768600193,
    .00000014340062895106,
    .00000000161384906966,
    .00000000001396650044,
    .00000000000009579451,
    .00000000000000053339,
    .00000000000000000245
  };
  static cheb_series bi0_cs = {
    bi0_data,
    11,
    -1, 1,
    11
  };

  static double ai0_data[21] = {
    .07575994494023796, 
    .00759138081082334,
    .00041531313389237,
    .00001070076463439,
    -.00000790117997921,
    -.00000078261435014,
    .00000027838499429,
    .00000000825247260,
    -.00000001204463945,
    .00000000155964859,
    .00000000022925563,
    -.00000000011916228,
    .00000000001757854,
    .00000000000112822,
    -.00000000000114684,
    .00000000000027155,
    -.00000000000002415,
    -.00000000000000608,
    .00000000000000314,
    -.00000000000000071,
    .00000000000000007
  };
  static cheb_series ai0_cs = {
    ai0_data,
    20,
    -1, 1,
    13
  };

  static double ai02_data[22] = {
    .05449041101410882,
    .00336911647825569,
    .00006889758346918,
    .00000289137052082,
    .00000020489185893,
    .00000002266668991,
    .00000000339623203,
    .00000000049406022,
    .00000000001188914,
    -.00000000003149915,
    -.00000000001321580,
    -.00000000000179419,
    .00000000000071801,
    .00000000000038529,
    .00000000000001539,
    -.00000000000004151,
    -.00000000000000954,
    .00000000000000382,
    .00000000000000176,
    -.00000000000000034,
    -.00000000000000027,
    .00000000000000003
  };
  static cheb_series ai02_cs = {
    ai02_data,
    21,
    -1, 1,
    11
  };

  double cheb_eval(const cheb_series * cs,
		   const double x)
  {
    int j;
    double d  = 0.0;
    double dd = 0.0;

    double y  = (2.0*x - cs->a - cs->b) / (cs->b - cs->a);
    double y2 = 2.0 * y;

    double e = 0.0;

    for(j = cs->order; j>=1; j--) {
      double temp = d;
      d = y2*d - dd + cs->c[j];
      e += fabs(y2*temp) + fabs(dd) + fabs(cs->c[j]);
      dd = temp;
    }

    { 
      double temp = d;
      d = y*d - dd + 0.5 * cs->c[0];
      e += fabs(y*temp) + fabs(dd) + 0.5 * fabs(cs->c[0]);
    }

    return d;
  }

  double gsl_sf_bessel_I0_scaled(const double x)
  {
    double y = fabs(x);

    if(y < 2.0 * GSL_SQRT_DBL_EPSILON) {
      return 1.0 - y;
    }
    else if(y <= 3.0) {
      const double ey = exp(-y);
      const double v = cheb_eval(&bi0_cs, y*y/4.5-1.0);
      return ey * (2.75 + v);
    }
    else if(y <= 8.0) {
      const double sy = sqrt(y);
      const double v = cheb_eval(&ai0_cs, (48.0/y-11.0)/5.0);
      return (0.375 + v) / sy;
    }
    else {
      const double sy = sqrt(y);
      const double v = cheb_eval(&ai02_cs, 16.0/y-1.0);
      return (0.375 + v) / sy;
    }
  }

  double gauss(double x,double sigma) {

    return std::exp(-x*x/(2*sigma*sigma));

   }

  double eval_gauss(double x, double sigma) {

    return 1./(2*M_PI*sigma*sigma)*std::exp(-(x*x)/(2*sigma*sigma));
  }
}

namespace Likelihood {

double SpatialGaussian::convolve(const Likelihood::MeanPsf& psf, 
				 double energy, double x, 
				 double sigma, int n) {

  const double xpmin = std::max(x - 10*sigma,0.0);
  const double xpmax = x + 10*sigma;
  const double dxp = (xpmax-xpmin)/double(n+1);
  const double s2 = sigma*sigma;
  
  double s = 0;
  for(int i = 0; i < n; i++)
    {
      double xp = xpmin + (double(i)+0.5)*dxp;
      double xx = x*xp/s2;
      double je = gsl_sf_bessel_I0_scaled(xx);      
      s += xp*psf(energy,xp)*je*std::exp(xx - (x*x+xp*xp)/(2*s2));
    }
  
  return s/s2*dxp;
}

SpatialGaussian::SpatialGaussian() 
  : SpatialFunction("SpatialGaussian",3) {
  m_width = 1.0;
  addParam("Width", 1.0, false);
  parameter("Width").setBounds(0.0, 180.);
}

SpatialGaussian::SpatialGaussian(double ra, double dec, double width) 
  : SpatialFunction("SpatialGaussian",3,ra,dec) {
  m_width = width;
  addParam("Width", m_width, false);
  parameter("Width").setBounds(0.0, 180.);
}

SpatialGaussian::SpatialGaussian(const SpatialGaussian & rhs) 
  : SpatialFunction(rhs), m_width(rhs.m_width) {

}

SpatialGaussian & SpatialGaussian::operator=(const SpatialGaussian & rhs) {
   if (this != &rhs) {
      SpatialFunction::operator=(rhs);
      m_width = rhs.m_width;
   }
   return *this;
}

SpatialGaussian::~SpatialGaussian() {
}

double SpatialGaussian::value(const astro::SkyDir & dir) const {
   double delta = this->dir().difference(dir)*180./M_PI;
   return gauss(delta,m_width)*std::pow(M_PI/180.,2);
}

double SpatialGaussian::value(double delta, double width) const {
   return gauss(delta,width)*std::pow(M_PI/180.,2);
}

double SpatialGaussian::value(const astro::SkyDir & dir, double energy, const MeanPsf& psf) const {
  double delta = dir.difference(this->dir())*180./M_PI;
  return SpatialGaussian::convolve(psf,energy,delta,m_width);
}

double SpatialGaussian::value(double delta, double energy, const MeanPsf& psf) const {
  return SpatialGaussian::convolve(psf,energy,delta,m_width);
}

void SpatialGaussian::update() {
  SpatialFunction::update();
  double width = getParam("Width").getValue();  
  m_width = width;
}

double SpatialGaussian::value(const optimizers::Arg & x) const {
   const SkyDirArg & dir = dynamic_cast<const SkyDirArg &>(x);
   double offset = dir().difference(this->dir())*180./M_PI;
   return value(offset,m_width);
}

double SpatialGaussian::derivByParamImp(const optimizers::Arg & x, 
                                      const std::string & parName) const {

   std::ostringstream message;
   message << "SpatialGaussian: cannot take derivative wrt "
           << "parameter " << parName;
   throw std::runtime_error(message.str());
}

} // namespace Likelihood
