/**
 * @file CountsMapHealpix.h
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/echarles/healpix_changes/Likelihood/Likelihood/CountsMapHealpix.h,v 1.4 2015/11/30 23:47:31 echarles Exp $
 */

#ifndef Likelihood_CountsMapHealpix_h
#define Likelihood_CountsMapHealpix_h

#include "Likelihood/CountsMapBase.h"

#include <string>

#include "astro/SkyDir.h"
#include "astro/HealpixProj.h"

#include "evtbin/DataProduct.h"
#include "evtbin/HealpixBinner.h"

#include "Likelihood/HistND.h"
#include "Likelihood/Pixel.h"

namespace astro {
   class SkyProj;
}

namespace evtbin {
   class Binner;
}

namespace tip {
   class Header;
}

namespace Likelihood {

/**
 * @class CountsMapHealpix
 * 
 */
   
class CountsMapHealpix : public CountsMapBase {

public:

   CountsMapHealpix(const std::string & countsMapFile);

   CountsMapHealpix(const CountsMapHealpix & counts_map);

   CountsMapHealpix(const CountsMapHealpix & counts_map,
		    unsigned int firstBin, unsigned int lastBin);

   virtual ~CountsMapHealpix() throw();

   virtual CountsMapBase* clone() const {
     return new CountsMapHealpix(*this);
   }

   virtual void binInput(tip::Table::ConstIterator begin, 
                         tip::Table::ConstIterator end);

   virtual void writeOutput(const std::string & creator, 
                            const std::string & out_file) const;
   
   virtual void setKeywords(tip::Header & header) const;

   virtual void getBoundaryPixelDirs(std::vector<astro::SkyDir> & pixelDirs) const;

   virtual bool withinBounds(const astro::SkyDir & dir, double energy, long border_size=0) const;

   virtual double pixelSize() const { return m_pixelSize; }

   inline int nPixels() const { return m_nPixels; }

   inline double solidAngle() const { return m_solidAngle; }

   inline bool allSky() const { return m_mapRadius >= 180.; }

   inline double mapRadius() const { return m_mapRadius; }

   inline const astro::HealpixProj* healpixProj() const { return m_healpixProj; }

   const std::vector<int>& pixelIndices() const;

   int globalToLocalIndex(int glo) const;
   
   int localToGlobalIndex(int loc) const;
      

protected:
   
   void latchCacheData();

private:

   // it is actually the same as m_proj
   // it is just here for convinience
   astro::HealpixProj* m_healpixProj;

   // it is actually the same as binners[0]
   // it is just here for convinience
   const evtbin::HealpixBinner* m_hpx_binner;

   double m_solidAngle;

   double m_pixelSize;
   
   double m_mapRadius;

   int m_nPixels;

   CountsMapHealpix & operator=(const CountsMapHealpix &) {return *this;}

   virtual double computeSolidAngle(std::vector<double>::const_iterator lon,
				    std::vector<double>::const_iterator lat,
				    const astro::ProjBase & proj) const { return m_solidAngle;}

   virtual void getPixels(std::vector<astro::SkyDir> & pixelDirs,
			  std::vector<double> & pixelSolidAngles) const;
   
   virtual void readKeywords(const std::string & countsMapFile);

   virtual void readImageData(const std::string & countsMapfile,
			      std::vector<evtbin::Binner *> & binners);

};

}

#endif // Likelihood_CountsMapHealpix_h
