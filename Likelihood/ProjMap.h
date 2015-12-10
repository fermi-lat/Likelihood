/**
 * @file ProjMap.h
 * @brief A map with reference point centered on the image and that
 * uses astro::ProjBase for indexing its internal representation.
 * @author E. Charles, from J. Chiang WcsMap2
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/echarles/healpix_changes/Likelihood/Likelihood/ProjMap.h,v 1.2 2015/03/03 05:59:56 echarles Exp $
 */

#ifndef Likelihood_ProjMap_h
#define Likelihood_ProjMap_h

#include <vector>
#include <utility>

#include "astro/SkyDir.h"
#include "astro/ProjBase.h"

namespace Likelihood {

class BinnedExposureBase;
class DiffuseSource;
class MeanPsf;

/**
 * @class ProjMap
 *
 */

class ProjMap {

public:

   ProjMap(const std::string & filename, 
           bool interpolate=true, bool enforceEnergyRange=false);

   virtual ~ProjMap();

   ProjMap(const ProjMap &);
 
   ProjMap(const ProjMap &, const double & energy);

   virtual ProjMap & operator=(const ProjMap &);

   virtual double operator()(const astro::SkyDir & dir, double energy=-1) const = 0;

   virtual double operator()(const astro::SkyDir & dir, int k) const = 0;

   virtual ProjMap* convolve(double energy, const MeanPsf & psf,
			     const BinnedExposureBase & exposure,
			     bool performConvolution=true,
			     int k=0) const = 0;
   
   virtual double solidAngle(double ilon, double ilat) const = 0;

   /// @return Pixel value as a function of index
   virtual double pixelValue(double ilon, double ilat, int k=0) const = 0;
   
   /// @return SkyDir corresponding to the pixel indices
   astro::SkyDir skyDir(double ilon, double ilat) const;  
   
   /// @return true if dir is within a cone half-angle m_mapRadius
   bool withinMapRadius(const astro::SkyDir & dir) const;

   int nenergies() const { return m_energies.size(); }

   const std::vector<double> & energies() const {
      return m_energies;
   }

   virtual std::pair<astro::SkyDir, astro::SkyDir> minMaxDistPixels(const astro::SkyDir & dir) const = 0;

   virtual bool insideMap(const astro::SkyDir & dir) const = 0;

   virtual double pixelSize() const = 0;
   
   inline double mapIntegral() const { return m_mapIntegral; }

   inline double mapRadius() const { return m_mapRadius; }

   double mapIntegral(double energy) const;

   /// Rebin the map data by the specified factor.  By default, this
   /// will average over combined pixels, as would be appropriate for
   /// intensity or exposure maps.  For counts maps, set
   /// average=false, so that combined bins are summed.  Return a new
   /// ProjMap object with the new geometry.  The reference direction
   /// will be unchanged from the original (and so will not generally
   /// point to the center of the map.)
   virtual ProjMap* rebin(unsigned int factor, bool average=true) = 0;

   inline void setExtrapolation(bool enforceEnergyRange) { m_enforceEnergyRange = enforceEnergyRange; }
   
   inline void setInterpolation(bool interpolate) { m_interpolate = interpolate; }

   inline bool enforceEnergyRange() const { return m_enforceEnergyRange; }

   inline unsigned long extrapolated() const { return m_extrapolated; }

   inline const astro::ProjBase* getProj() const { return m_proj; }

   static double interpolatePowerLaw(double x, double x1, double x2,
                                     double y1, double y2);

protected:

   ProjMap();

   // Query private data
   inline const astro::SkyDir& getRefDir() const { return m_refDir; }
   
   inline bool getInterpolate() const { return m_interpolate; }

   // Set and access private data
   void setProjInfo(const astro::SkyDir& dir, const astro::ProjBase& proj);

   inline void setMapRadius(const double& val) { m_mapRadius = val; }

   inline void setFilename(const std::string& filename) { 
     m_filename = filename;
   }

   inline std::vector<double> & energies_access() {
     return m_energies;
   }

   inline double& mapIntegral_access() { return m_mapIntegral; }

   inline std::vector<float>& mapIntegrals_access() { return m_mapIntegrals; }

   inline unsigned long& extrapolated_access() const { return m_extrapolated; }
       
   virtual void computeMapIntegrals() = 0;

   void check_energy_index(int k) const;

   void check_energy(double energy) const;

private:

   std::string m_filename;
 
   astro::SkyDir m_refDir;

   // enclosing cone half angle for map
   double m_mapRadius;

   astro::ProjBase* m_proj;
   
   std::vector<double> m_energies;

   bool m_interpolate;

   bool m_enforceEnergyRange;
 
   mutable unsigned long m_extrapolated;

   double m_mapIntegral;

   std::vector<float> m_mapIntegrals;


};

} // namespace Likelihood

#endif // Likelihood_ProjMap_h
