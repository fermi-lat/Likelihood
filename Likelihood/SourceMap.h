/**
 * @file SourceMap.h
 * @brief Spatial distribution of a source folded through the
 *        instrument response.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceMap.h,v 1.6 2004/10/07 00:01:02 jchiang Exp $
 */

#ifndef Likelihood_SourceMap_h
#define Likelihood_SourceMap_h

#include "Likelihood/BinnedExposure.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/Pixel.h"

namespace Likelihood {

class Source;
class CountsMap;

/*
 * @class SourceMap
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceMap.h,v 1.6 2004/10/07 00:01:02 jchiang Exp $
 */

class SourceMap {

public:

   SourceMap(Source * src, const CountsMap & dataMap);

   SourceMap(const std::string & sourceMapsFile, const std::string & srcName);

   ~SourceMap() {
      delete s_meanPsf;
      delete s_binnedExposure;
   }

   const std::vector<double> & model() const {return m_model;}

   const std::vector<double> & npreds() const {return m_npreds;}
   
   static void setBinnedExposure(const std::string & filename) {
      s_binnedExposure = new BinnedExposure(filename);
   }

private:

   std::string m_name;

   const CountsMap & m_dataMap;

/// @brief m_models has the same size as the data in the dataMap plus
/// one energy plane.
///
/// @todo Keep track of event types included in a given SourceMap.
   std::vector<double> m_model;

/// @brief Each entry is the angular integral over the energy plane.
   std::vector<double> m_npreds;

   class Aeff : public Pixel::Aeff {
   public:
      Aeff(Source * src, const astro::SkyDir & appDir,
           double energy, int type)
         : Pixel::Aeff(src, appDir, energy, type) {}
      virtual double operator()(double costheta) const;
   };

   static MeanPsf * s_meanPsf;
   static BinnedExposure * s_binnedExposure;

   static std::vector<double> s_phi;
   static std::vector<double> s_mu;

   double sourceRegionIntegral(Source * src, const Pixel & pixel,
                               double energy) const;

   void prepareAngleArrays(int nmu, int nphi);

   void getCelestialDir(double phi, double mu, 
                        FitsImage::EquinoxRotation & eqRot,
                        astro::SkyDir & dir) const;
};

} // namespace Likelihood

#endif // Likelihood_SourceMap_h
