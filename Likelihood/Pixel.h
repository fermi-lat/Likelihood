/**
 * @file Pixel.h
 * @brief Provide access to predicted numbers of counts in a pixel and
 * derivatives wrt model parameters.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef Likelihood_Pixel_h
#define Likelihood_Pixel_h

#include <string>
#include <vector>

#include "astro/SkyDir.h"
#include "map_tools/Exposure.h"
//#include "Likelihood/SourceModel.h"

namespace Likelihood {

   class Source;
   class SourceModel;

/**
 * @class Pixel
 *
 * @brief A Flyweight-like class that provides access to model
 * predictions for numbers of counts and derivatives of the counts
 * prediction wrt model parameters.  The extrinsic quantities passed
 * to the public methods are the bounds of the energy band of interest
 * (see GOF, p. 195)
 *
 * @author J. Chiang
 *
 * $Header$
 */

class Pixel {

public:

   Pixel(double ra=0, double dec=0, double solidAngle=0) 
      : m_dir(astro::SkyDir(ra, dec)), m_solidAngle(solidAngle) {}

   Pixel(const astro::SkyDir & dir, double solidAngle) 
      : m_dir(dir), m_solidAngle(solidAngle) {}

   double modelCounts(double emin, double emax) const;

   void getFreeDerivs(double emin, double emax, 
                      std::vector<double> & derivs) const;

   class Aeff : public map_tools::Exposure::Aeff {
   public:
      Aeff(Source * src, const astro::SkyDir & appDir, 
           double energy, int type);
      virtual double operator()(double costheta) const;
   private:
      Source * m_src;
      const astro::SkyDir & m_appDir;
      double m_energy;
      int m_type;
      double m_separation;
   };

   class AeffDeriv : public map_tools::Exposure::Aeff {
   public:
      AeffDeriv(Source * src, const std::string & paramName, 
                const astro::SkyDir & appDir, double energy, int type);
      virtual ~AeffDeriv() {}
      virtual double operator()(double costheta) const;
   private:
      Source * m_src;
      std::string m_paramName;
      const astro::SkyDir & m_appDir;
      double m_energy;
      int m_type;
      double m_separation;
   };

private:

   astro::SkyDir m_dir;
   double m_solidAngle;

   static SourceModel * s_srcModel;
   
};

} // namespace Likelihood

#endif // Likelihood_Pixel_h
