/**
 * @file MapBase.h
 * @brief Base class for FITS map objects.
 * 
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/MapBase.h,v 1.1 2009/02/18 06:57:42 jchiang Exp $
 */

#ifndef Likelihood_MapBase_h
#define Likelihood_MapBase_h

#include <string>
#include <utility>
#include <vector>

#include "astro/SkyDir.h"

namespace Likelihood {

class WcsMap;

/**
 * @class MapBase
 * @brief Base class for FITS map objects.
 */

class MapBase {

public:

   MapBase();
   
   MapBase(const std::string & fitsFile, const std::string & extension="");

   MapBase(const MapBase & other);

   MapBase & operator=(const MapBase & rhs);

   virtual ~MapBase();

   virtual void readFitsFile(const std::string & fitsFile,
                             const std::string & extension="");

   virtual bool insideMap(const astro::SkyDir & dir) const;

   virtual void getDiffRespLimits(const astro::SkyDir & dir, 
                                  double & mumin, double & mumax,
                                  double & phimin, double & phimax) const;

protected:

   WcsMap * m_wcsmap;

   std::string m_fitsFile;
   std::string m_extension;

private:
   
   void getMinMaxDistPixels(const astro::SkyDir &,
                            astro::SkyDir & closestPixel, 
                            astro::SkyDir & farthestPixel) const;
   
   void getCorners(std::vector<astro::SkyDir> & corners) const;

};

} // namespace Likelihood

#endif // Likelihood_MapBase_h
