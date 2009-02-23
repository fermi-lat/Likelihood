/**
 * @file MapBase.h
 * @brief Base class for FITS map objects.
 * 
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/MapBase.h,v 1.3 2009/02/21 02:03:18 jchiang Exp $
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

   virtual double mapIntegral(double energy) const = 0;

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

class MapBaseException : public std::exception {
public:
   MapBaseException() {}
   MapBaseException(std::string errorString) : m_what(errorString) {}
   virtual ~MapBaseException() throw() {}
   virtual const char *what() const throw() {return m_what.c_str();}
protected:
   std::string m_what;
};

} // namespace Likelihood

#endif // Likelihood_MapBase_h
