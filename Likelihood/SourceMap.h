/**
 * @file SourceMap.h
 * @brief Spatial distribution of a source folded through the
 *        instrument response.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceMap.h,v 1.1 2004/09/22 05:39:57 jchiang Exp $
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceMap.h,v 1.1 2004/09/22 05:39:57 jchiang Exp $
 */

class SourceMap {

public:

   SourceMap(Source * src, const CountsMap & dataMap);

   ~SourceMap() {}

   const std::vector<double> & model() const {return m_model;}

private:

/// m_models has the same size as the data in the dataMap plus one
/// energy plane.
///
/// @todo Keep track of event types included in a given SourceMap.

   std::vector<double> m_model;

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
