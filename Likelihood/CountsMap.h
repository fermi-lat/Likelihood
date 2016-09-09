/**
 * @file CountsMap.h
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/CountsMap.h,v 1.33 2015/12/10 00:57:58 echarles Exp $
 */

#ifndef Likelihood_CountsMap_h
#define Likelihood_CountsMap_h

// EAC, make a base class for CountsMap
#include "Likelihood/CountsMapBase.h"

#include <string>

#include "astro/SkyDir.h"

#include "evtbin/DataProduct.h"

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
 * @class CountsMap
 * 
 */
   
class CountsMap : public CountsMapBase {

public:

   CountsMap(const std::string & event_file, const std::string & ev_table,
             const std::string & sc_file, const std::string & sc_table,
             double ref_ra, double ref_dec, const std::string & proj,
             unsigned long num_x_pix, unsigned long num_y_pix, 
             double pix_scale, double axis_rot, bool use_lb, 
             const std::string & ra_field, const std::string & dec_field,
             double emin, double emax, unsigned long nenergies);

   CountsMap(const std::string & event_file, const std::string & ev_table,
             const std::string & sc_file, const std::string & sc_table,
             double ref_ra, double ref_dec, const std::string & proj,
             unsigned long num_x_pix, unsigned long num_y_pix, 
             double pix_scale, double axis_rot, bool use_lb, 
             const std::string & ra_field, const std::string & dec_field,
             const std::vector<double> & energies);

   CountsMap(const std::string & event_file, const std::string & ev_table,
             const std::string & sc_file, const std::string & sc_table,
             double ref_ra, double ref_dec, const std::string & proj,
             unsigned long num_x_pix, unsigned long num_y_pix, 
             double x_pix_scale, double y_pix_scale, 
             double axis_rot, bool use_lb, 
             const std::string & ra_field, const std::string & dec_field,
             const std::vector<double> & energies);

   CountsMap(const std::string & event_file, const std::string & ev_table,
             const std::string & sc_file, const std::string & sc_table,
             double ref_ra, double ref_dec, const std::string & proj,
             unsigned long num_x_pix, unsigned long num_y_pix, 
             double x_pix_scale, double y_pix_scale, 
             double axis_rot, bool use_lb, 
             const std::string & ra_field, const std::string & dec_field,
             const std::vector<double> & emins, 
             const std::vector<double> & emaxs);

   CountsMap(const std::string & countsMapFile);

   CountsMap(const CountsMap & counts_map);

   CountsMap(const CountsMap & counts_map, unsigned int firstBin, unsigned int lastBin);

   virtual ~CountsMap() throw();

   virtual CountsMapBase* clone() const { 
     return new CountsMap(*this);
   }

   virtual void binInput(tip::Table::ConstIterator begin, 
                         tip::Table::ConstIterator end);

   virtual void writeOutput(const std::string & creator, 
                            const std::string & out_file) const;
   
   virtual void setKeywords(tip::Header & header) const;

   virtual void getBoundaryPixelDirs(std::vector<astro::SkyDir> & pixelDirs) const;

   virtual bool withinBounds(const astro::SkyDir & dir, double energy, long border_size=0) const;

   virtual double pixelSize() const;

   virtual double mapRadius() const { return 0.7071067811865476 * std::max( -1. * m_naxes[0] * m_cdelt[0],  m_naxes[1] * m_cdelt[1] ); }

   double cdelt1() const {return m_cdelt[0];}
   double cdelt2() const {return m_cdelt[1];}
   double crpix1() const {return m_crpix[0];}
   double crpix2() const {return m_crpix[1];}
   double crval1() const {return m_crval[0];}
   double crval2() const {return m_crval[1];}
   double crota2() const {return m_axis_rot;}
   long naxis1() const {return m_naxes[0];}
   long naxis2() const {return m_naxes[1];}

   bool conformingMap() const {return m_conforms;}

protected:

   long m_naxes[3];
   double m_crpix[3];
   double m_crval[3];
   double m_cdelt[3];
   double m_axis_rot;

private:

   // Flag to indicate that this counts map conforms to those created by
   // gtbin, i.e., with the reference pixel in the center and 
   // with cdelt1==cdelt2.
   bool m_conforms;

   CountsMap & operator=(const CountsMap &) {return *this;}

   void init(std::vector<evtbin::Binner *> & binners, 
             const std::string & event_file, const std::string & ev_table,
             const std::string & sc_file, const std::string & sc_table,
             unsigned long num_x_pix, unsigned long num_y_pix, 
             double ref_ra, double ref_dec, double x_pix_scale, 
             double y_pix_scale, double emin,
             double emax, unsigned long nenergies, bool use_lb, 
             const std::string & proj);

   virtual double computeSolidAngle(std::vector<double>::const_iterator lon,
				    std::vector<double>::const_iterator lat,
				    const astro::ProjBase & proj) const;

   virtual void getPixels(std::vector<astro::SkyDir> & pixelDirs,
			  std::vector<double> & pixelSolidAngles) const;
   
   virtual void readKeywords(const std::string & countsMapFile);

   virtual void readImageData(const std::string & countsMapfile,
			      std::vector<evtbin::Binner *> & binners);

   void checkMapConforms();
};

}

#endif // Likelihood_CountsMap_h
