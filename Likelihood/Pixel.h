/**
 * @file Pixel.h
 * @brief Provide access to predicted numbers of counts in a pixel and
 * derivatives wrt model parameters.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/Pixel.h,v 1.12 2012/03/28 22:00:42 jchiang Exp $
 */

#ifndef Likelihood_Pixel_h
#define Likelihood_Pixel_h

#include <string>
#include <vector>

#include "astro/SkyDir.h"

#include "Likelihood/ExposureCube.h"

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
 * and the SourceModel object (GOF, p. 195).
 *
 */

class Pixel {

public:

   Pixel(double ra=0, double dec=0, double solidAngle=0) 
      : m_dir(astro::SkyDir(ra, dec)), m_solidAngle(solidAngle) {}

   Pixel(const astro::SkyDir & dir, double solidAngle) 
      : m_dir(dir), m_solidAngle(solidAngle) {}

   double modelCounts(double emin, double emax, SourceModel & srcModel) const;

   void getFreeDerivs(double emin, double emax, SourceModel & srcModel,
                      std::vector<double> & derivs) const;

   const astro::SkyDir & dir() const {return m_dir;}

   double solidAngle() const {return m_solidAngle;}

   /// @return The iterator to the Pixel object that is closest to the
   /// input Pixel object in angular distance.
   template<class InIt>
   static InIt find(InIt first, InIt last, const Pixel & val, 
                    double fudgeFactor=1) {
      InIt closest = first;
      double dist(first->dir().difference(val.dir()));
      for (InIt candidate = first; candidate != last; ++candidate) {
         double newdist(candidate->dir().difference(val.dir()));
         if (newdist < dist) {
            dist = newdist;
            closest = candidate;
         }
      }
/// @todo Use wcslib information to determine Pixel shape and whether
/// the current point is really inside the nearest Pixel. Here we use
/// a crude approximation.
      if ((1. - closest->solidAngle()/2./M_PI) < std::cos(dist/fudgeFactor)) {
         return closest;
      }
      return last;
   }

#ifndef SWIG
   class Aeff : public ExposureCube::AeffBase {
   public:
      Aeff(Source * src, const astro::SkyDir & appDir, 
           double energy, int type, double time);
   protected:
      Source * m_src;
      const astro::SkyDir & m_appDir;
      double m_energy;
      int m_type;
      double m_separation;
      double m_time;

      virtual double value(double costheta, double phi=0) const;
   };

   class AeffDeriv : public ExposureCube::AeffBase {
   public:
      AeffDeriv(Source * src, const std::string & paramName, 
                const astro::SkyDir & appDir, double energy, int type,
                double time);
      virtual ~AeffDeriv() {}
   protected:
      Source * m_src;
      std::string m_paramName;
      const astro::SkyDir & m_appDir;
      double m_energy;
      int m_type;
      double m_separation;
      double m_time;

      virtual double value(double costheta, double phi=0) const;
   };
#endif // SWIG

private:

   astro::SkyDir m_dir;
   double m_solidAngle;

};

} // namespace Likelihood

#endif // Likelihood_Pixel_h
