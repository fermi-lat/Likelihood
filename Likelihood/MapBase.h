/**
 * @file MapBase.h
 * @brief Base class for FITS map objects.
 * 
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/MapBase.h,v 1.8 2011/03/15 05:37:31 jchiang Exp $
 */

#ifndef Likelihood_MapBase_h
#define Likelihood_MapBase_h

#include <string>
#include <utility>
#include <vector>

#include "astro/SkyDir.h"

#include "Likelihood/WcsMap2.h"

namespace Likelihood {

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
                             const std::string & extension="",
                             bool loadMap=true);

   virtual void readFitsFile();

   virtual bool insideMap(const astro::SkyDir & dir) const;

   virtual void getDiffRespLimits(const astro::SkyDir & dir, 
                                  double & mumin, double & mumax,
                                  double & phimin, double & phimax) const;

   virtual const std::string & fitsFile() const {
      return m_fitsFile;
   }

   const WcsMap2 & wcsmap() const {
      if (!m_wcsmap) {
         const_cast<MapBase *>(this)->readFitsFile();
      }
      return *m_wcsmap;
   }

   virtual WcsMap2 & wcsmap();

   virtual void deleteMap() {
      delete m_wcsmap;
      m_wcsmap = 0;
   }

   virtual void rebin(unsigned int factor, bool average=true);

protected:

   WcsMap2 * m_wcsmap;

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
