/**
 * @file CountsMap.cxx
 *
 * $Header$
 */

#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "astro/SkyDir.h"
#include "astro/SkyProj.h"

#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Table.h"
#include "tip/tip_types.h"

#include "evtbin/LinearBinner.h"
#include "evtbin/LogBinner.h"

#include "Likelihood/CountsMap.h"
#include "Likelihood/HistND.h"

namespace Likelihood {

CountsMap::CountsMap(const std::string & event_file, 
                     const std::string & sc_file, 
                     double ref_ra, double ref_dec, const std::string & proj, 
                     unsigned long num_x_pix, unsigned long num_y_pix, 
                     double pix_scale, double axis_rot, bool use_lb,
                     const std::string & ra_field, 
                     const std::string & dec_field, 
                     double emin, double emax, unsigned long nenergies) 
   : DataProduct(event_file), m_hist(0), m_proj_name(proj), 
     m_crpix(), m_crval(), m_cdelt(), m_axis_rot(axis_rot), 
     m_use_lb(use_lb), m_proj(0) {

// These binners will be deleted by the Hist base class.
   std::vector<Binner *> binners;

// The astro::SkyDir::project method will convert ra and dec into 
// bin indices, so we set up the LinearBinners generically by index.
   binners.push_back(new LinearBinner(0.5, num_x_pix + 0.5, 1., ra_field));
   binners.push_back(new LinearBinner(0.5, num_y_pix + 0.5, 1., dec_field));

// Energy bins are set up normally.
   binners.push_back(new LogBinner(emin, emax, nenergies, "photon energy"));

   m_hist = new HistND(binners);

   m_naxes[0] = num_x_pix;
   m_naxes[1] = num_y_pix;
   m_naxes[2] = nenergies;
   m_crpix[0] = (num_x_pix + 1.) / 2.;
   m_crpix[1] = (num_y_pix + 1.) / 2.;
   m_crpix[2] = 1;
   m_crval[0] = ref_ra;
   m_crval[1] = ref_dec;
   m_crval[2] = emin;
   m_cdelt[0] = -pix_scale;
   m_cdelt[1] = pix_scale;
   m_cdelt[2] = log(emax/emin)/(nenergies-1.);
   
   m_proj = new astro::SkyProj(proj, m_crpix, m_crval, m_cdelt, 
                               m_axis_rot, use_lb);

// Collect any/all needed keywords from the event file.
   harvestKeywords(event_file, "EVENTS");

// Correct time keywords.
   adjustTimeKeywords(sc_file);
}

CountsMap::CountsMap(const CountsMap & rhs) {
   m_hist = rhs.m_hist->clone();
   m_proj_name = rhs.m_proj_name;
   for (int i = 0; i < 3; i++) {
      m_naxes[i] = rhs.m_naxes[i];
      m_crpix[i] = rhs.m_crpix[i];
      m_crval[i] = rhs.m_crval[i];
      m_cdelt[i] = rhs.m_cdelt[i];
   }
   m_axis_rot = rhs.m_axis_rot;
   m_use_lb = rhs.m_use_lb;
   m_proj = new astro::SkyProj(m_proj_name, m_crpix, m_crval, m_cdelt, 
                               m_axis_rot, m_use_lb);
}

CountsMap::~CountsMap() throw() { 
   try {
      delete m_proj;
      delete m_hist;
   } catch (std::exception & eObj) {
      std::cerr << eObj.what() << std::endl;
   } catch (...) {
   }
}

void CountsMap::binInput(tip::Table::ConstIterator begin, 
                         tip::Table::ConstIterator end) {

// Get binners for the three dimensions.
   const Hist::BinnerCont_t & binners = m_hist->getBinners();

// From each sky binner, get the name of its field, interpreted as ra
// and dec.
   std::string ra_field = binners[0]->getName();
   std::string dec_field = binners[1]->getName();
   
// Fill histogram, converting each RA/DEC to Sky X/Y on the fly:

   for (tip::Table::ConstIterator itor = begin; itor != end; ++itor) {
      // Extract the ra and dec from each record.
      double ra = (*itor)[ra_field].get();
      double dec = (*itor)[dec_field].get();
      
      // Convert to sky coordinates.
      std::pair<double, double> coord 
         = astro::SkyDir(ra, dec).project(*m_proj);

      // Get the event energy.
      double energy = (*itor)["ENERGY"].get();

      double my_values[] = {coord.first, coord.second, energy};
      std::vector<double> values(my_values, my_values + 3);

      // Bin the value.
      m_hist->fillBin(values);
   }
}

void CountsMap::writeOutput(const std::string & creator, 
                            const std::string & out_file) const {

// Standard file creation from base class.
    createFile(creator, out_file, m_data_dir + "LatCountsMapTemplate");

// Open Count map extension of output PHA1 file. Use an auto_ptr so
// that the table object will for sure be deleted, even if an
// exception is thrown.
    std::auto_ptr<tip::Image> 
       output_image(tip::IFileSvc::instance().editImage(out_file, ""));

// Get dimensions of image.
    typedef std::vector<tip::PixOrd_t> DimCont_t;
    DimCont_t dims = output_image->getImageDimensions();

// Make sure image is three dimensional.
    DimCont_t::size_type num_dims = dims.size();
    if (3 != num_dims) {
       throw std::runtime_error("CountsMap::writeOutput "
                                + std::string("cannot write a count map ")
                                + "to an image which is not 3D");
    }

// Get the binners.
    const Hist::BinnerCont_t & binners = m_hist->getBinners();

// Write c* keywords
    tip::Header & header = output_image->getHeader();
    header["CRPIX1"].set(m_crpix[0]);
    header["CRPIX2"].set(m_crpix[1]);
    header["CRPIX3"].set(m_crpix[2]);
    header["CRVAL1"].set(m_crval[0]);
    header["CRVAL2"].set(m_crval[1]);
    header["CRVAL3"].set(m_crval[2]);
    header["CDELT1"].set(m_cdelt[0]);
    header["CDELT2"].set(m_cdelt[1]);
    header["CDELT3"].set(m_cdelt[2]);
    header["CROTA2"].set(m_axis_rot);
    header["CTYPE1"].set(binners[0]->getName() + "---" + m_proj_name);
    header["CTYPE2"].set(binners[1]->getName() + "---" + m_proj_name);
    header["CTYPE3"].set(binners[2]->getName());

// Resize image dimensions to conform to the binner dimensions.
    for (DimCont_t::size_type index = 0; index != num_dims; ++index) {
       dims[index] = binners.at(index)->getNumBins();
    }

// Set size of image.
    output_image->setImageDimensions(dims);

// Copy bins into image.
    std::vector<float> float_image(m_hist->data().size());
    std::copy(m_hist->data().begin(), m_hist->data().end(),
              float_image.begin());
    output_image->set(float_image);

// Write the GTI extension.
    writeGti(out_file);
//    std::cout << "Done." << std::endl;
}

void CountsMap::setImage(const std::vector<double> & image) {
   m_hist->setData(image);
}


long CountsMap::imageDimension(int idim) const {
   const Hist::BinnerCont_t & binners = m_hist->getBinners();
   if (i < 0 || i > 2) {
      throw std::invalid_argument(std::string("CountsMap::imageDimension:\n")
                                  + "Invalid image dimension value.");
   }
   return binners[i]->getNumBins();
}

void CountsMap::getAxisVector(int idim, std::vector<double> axisVector) const {
   const Hist::BinnerCont_t & binners = m_hist->getBinners();
   if (i < 0 || i >2) {
      throw std::invalid_argument(std::string("CountsMap::getAxisVector:\n")
                                  + "Invalid image dimension value.");
   }
   axisVector.clear();
   for (long j = 0; j < binners[i]->getNumBins(); j++) {
      axisVector.push_back(binners[i]->getInterval(j).begin());
   }
   long j = binners[i]->getNumBins() - 1;
   axisVector.push_back(binners[i]->getInterval(j).end());
}

} // namespace Likelihood
