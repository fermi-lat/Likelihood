/**
 * @file SourceMap.h
 * @brief Spatial distribution of a source folded through the
 *        instrument response.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef Likelihood_SourceMap_h
#define Likelihood_SourceMap_h

#include "Likelihood/Pixel.h"

namespace Likelihood {

class Source;
class CountsMap;

/*
 * @class SourceMap
 *
 * $Header$
 */

class SourceMap {

public:

   SourceMap(Source * src, const CountsMap & dataMap);

   ~SourceMap() {}

   const std::vector<double> & model(int evtType) const {
      return m_models.at(evtType);
   }

private:

// Each vector in m_models has the same size as the data in the
// dataMap plus one energy plane.  These vectors of vectors are
// indexed by event types.
//
// @todo Ensure this implementation can accommodate the possible
// representations of event types
   std::vector< std::vector<double> > m_models;

   std::string m_name;

   class Aeff : public Pixel::Aeff {
   public:
      Aeff(Source * src, const astro::SkyDir & appDir,
           double energy, int type)
         : Pixel::Aeff(src, appDir, energy, type) {}
      virtual double operator()(double costheta) const;
   };
};

} // namespace Likelihood

#endif // Likelihood_SourceMap_h
