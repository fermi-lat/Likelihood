/**
 * @file CountsMap.h
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/CountsMap.h,v 1.4 2004/09/03 06:08:56 jchiang Exp $
 */

#ifndef Likelihood_CountsMap_h
#define Likelihood_CountsMap_h

#include <string>

#include "evtbin/DataProduct.h"

#include "Likelihood/HistND.h"

namespace astro {
   class SkyProj;
}

namespace Likelihood {

#include "fitsio.h"

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

};

}

#endif // Likelihood_CountsMap_h
