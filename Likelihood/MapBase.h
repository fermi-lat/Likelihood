/**
 * @file MapBase.h
 * @brief Base class for FITS map objects.
 * 
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/echarles/healpix_changes/Likelihood/Likelihood/MapBase.h,v 1.4 2015/03/05 19:58:25 echarles Exp $
 */

#ifndef Likelihood_MapBase_h
#define Likelihood_MapBase_h

#include <string>
#include <utility>
#include <vector>

#include "astro/SkyDir.h"

#include "Likelihood/ProjMap.h"

namespace Likelihood {

class ExposureMap;   
// EAC, add projection specific methods 
class WcsMap2;
class HealpixProjMap;

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

   void getDiffRespLimits_wcs(const astro::SkyDir & dir,
			      const WcsMap2& wcsmap,
			      double & mumin, double & mumax,
			      double & phimin, double & phimax) const;
   
   void getDiffRespLimits_healpix(const astro::SkyDir & dir,
				  const HealpixProjMap& healmap,
				  double & mumin, double & mumax,
				  double & phimin, double & phimax) const;


   virtual const std::string & fitsFile() const {
      return m_fitsFile;
   }

   // EAC, switch to using ProjMap base class
   const ProjMap & projmap() const {
      if (!m_projmap) {
         const_cast<MapBase *>(this)->readFitsFile();
      }
      return *m_projmap;
   }

   // EAC, switch to using ProjMap base class
   virtual ProjMap & projmap();

   virtual void deleteMap();

   // Update the status of the contained WcsMap2 pointer 
   // (typically in response to WcsMapLibrary::notify())
   virtual void update();

   virtual void rebin(unsigned int factor, bool average=true);

   virtual void integrateSpatialDist(const std::vector<double> & energies,
                                     const ExposureMap & expmap,
                                     std::vector<double> & exposure) const = 0;

protected:

   std::string m_fitsFile;
   std::string m_expandedFileName;
   std::string m_extension;

   static double interpolatePowerLaw(double x, double x1, double x2,
                                     double y1, double y2);

private:

   // EAC, switch to using ProjMap base class
   ProjMap * m_projmap;
   
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
