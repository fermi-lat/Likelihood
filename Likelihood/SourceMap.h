/**
 * @file SourceMap.h
 * @brief Spatial distribution of a source folded through the
 *        instrument response.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceMap.h,v 1.20 2005/02/27 06:42:24 jchiang Exp $
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceMap.h,v 1.20 2005/02/27 06:42:24 jchiang Exp $
 */

class SourceMap {

public:

   SourceMap(Source * src, const CountsMap * dataMap,
             const Observation & observation);

   SourceMap(const std::string & sourceMapsFile, const std::string & srcName);

   ~SourceMap() {
      s_refCount--;
      if (s_refCount == 0) {
         delete s_meanPsf;
         s_meanPsf = 0;
         delete s_binnedExposure;
         s_binnedExposure = 0;
      }
      if (m_deleteDataMap) delete m_dataMap;
   }

   const std::vector<double> & model() const {return m_model;}

   const std::vector<double> & npreds() const {return m_npreds;}
   
   static void setBinnedExposure(const std::string & filename) {
      s_binnedExposure = new BinnedExposure(filename);
   }

   void save(const std::string & filename) const;

private:

   std::string m_name;

   const CountsMap * m_dataMap;

   bool m_deleteDataMap;

/// @brief m_models has the same size as the data in the dataMap plus
/// one energy plane.
///
/// @todo Keep track of event types included in a given SourceMap.
   std::vector<double> m_model;

/// @brief Each entry is the angular integral over the energy plane.
   std::vector<double> m_npreds;

   std::vector<double> m_energies;

/// @brief This vector of SkyDir objects is used by
/// sourceRegionIntegral for diffuse sources
   std::vector<astro::SkyDir> m_srcDirs;

   std::vector<double> m_srcStrengths;

   class Aeff : public Pixel::Aeff {
   public:
      Aeff(Source * src, const astro::SkyDir & appDir,
           double energy, int type, const Observation & observation)
         : Pixel::Aeff(src, appDir, energy, type),
           m_observation(observation) {}
      virtual double operator()(double costheta) const;
   private:
      const Observation & m_observation;
   };

   static MeanPsf * s_meanPsf;
   static BinnedExposure * s_binnedExposure;
   static unsigned int s_refCount;

   static std::vector<double> s_phi;
   static std::vector<double> s_mu;
   static std::vector<double> s_theta;

   double sourceRegionIntegral(double energy,
                               const Observation & observation) const;

   void computeSrcDirs(const Pixel & pixel, Source * src);

   void prepareAngleArrays(int nmu=100, int nphi=50);

   void getCelestialDir(double phi, double mu, 
                        FitsImage::EquinoxRotation & eqRot,
                        astro::SkyDir & dir) const;

   void fitsReportError(FILE *stream, int status) const;

   bool haveMapCubeFunction(DiffuseSource * src) const;

   void recomputeSrcStrengths(DiffuseSource * src, double energy);

};

} // namespace Likelihood

#endif // Likelihood_SourceMap_h
