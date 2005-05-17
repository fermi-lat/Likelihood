/**
 * @file CountsMap.h
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/CountsMap.h,v 1.14 2005/05/14 00:59:09 jchiang Exp $
 */

#ifndef Likelihood_CountsMap_h
#define Likelihood_CountsMap_h

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
   
class CountsMap : public evtbin::DataProduct {

public:

   CountsMap(const std::string & event_file, const std::string & sc_file,
             double ref_ra, double ref_dec, const std::string & proj,
             unsigned long num_x_pix, unsigned long num_y_pix, 
             double pix_scale, double axis_rot, bool use_lb, 
             const std::string & ra_field, const std::string & dec_field,
             double emin, double emax, unsigned long nenergies);

   CountsMap(const std::string & event_file, const std::string & sc_file,
             double ref_ra, double ref_dec, const std::string & proj,
             unsigned long num_x_pix, unsigned long num_y_pix, 
             double pix_scale, double axis_rot, bool use_lb, 
             const std::string & ra_field, const std::string & dec_field,
             const std::vector<double> & energies);

   CountsMap(const std::string & countsMapFile);

   CountsMap(const CountsMap & counts_map);

   virtual ~CountsMap() throw();

   virtual void binInput(tip::Table::ConstIterator begin, 
                         tip::Table::ConstIterator end);

   virtual void writeOutput(const std::string & creator, 
                            const std::string & out_file) const;

   void setImage(const std::vector<double> & image);
   
   long imageDimension(int idim) const;

   void getAxisVector(int idim, std::vector<double> & axisVector) const;

   const astro::SkyProj & projection() const {return *m_proj;}

   const std::vector<double> & data() const {
      return m_hist->data();
   }

   void setKeywords(tip::Header & header) const;

   void getPixels(std::vector<Pixel> & pixels) const;

   void getBoundaryPixelDirs(std::vector<astro::SkyDir> & pixelDirs) const;

   const astro::SkyDir & mapCenter() const {return m_center;}

   bool withinBounds(const astro::SkyDir & dir, double energy) const;

protected:

   HistND * m_hist;
   std::string m_proj_name;
   long m_naxes[3];
   double m_crpix[3];
   double m_crval[3];
   double m_cdelt[3];
   double m_axis_rot;
   bool m_use_lb;
   astro::SkyProj * m_proj;

private:

   astro::SkyDir m_center;

   CountsMap & operator=(const CountsMap & rhs) {return *this;}

   void init(std::vector<evtbin::Binner *> & binners, 
             const std::string & event_file, 
             const std::string & sc_file, unsigned long num_x_pix, 
             unsigned long num_y_pix, 
             double ref_ra, double ref_dec, double pix_scale, double emin,
             double emax, unsigned long nenergies, bool use_lb, 
             const std::string & proj);

   double computeSolidAngle(std::vector<double>::const_iterator lon,
                            std::vector<double>::const_iterator lat,
                            const astro::SkyProj & proj) const;

   void getPixels(std::vector<astro::SkyDir> & pixelDirs,
                  std::vector<double> & pixelSolidAngles) const;

   void readKeywords(const std::string & countsMapFile);
   void readEbounds(const std::string & countsMapfile,
                    std::vector<evtbin::Binner *> & binners);
   void readImageData(const std::string & countsMapfile,
                      std::vector<evtbin::Binner *> & binners);

   void setCenter();

   void setDataDir();

   void deleteBinners(std::vector<evtbin::Binner *> & binners) const;
};

}

#endif // Likelihood_CountsMap_h
