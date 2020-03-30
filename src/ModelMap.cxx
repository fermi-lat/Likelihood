/**
 * @file ModelMap.cxx
 * @brief Standalone class for compute a model map using a 
 * BinnedLikelihood object.
 */

#include <memory>
#include <sstream>
#include <stdexcept>
#include <ctime>

#include "fitsio.h"

#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Table.h"

#include "st_facilities/FitsUtil.h"

#include "dataSubselector/Cuts.h"

#include "Likelihood/BinnedLikelihood.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/CountsMapHealpix.h"

#include "Likelihood/ModelMap.h"

namespace {
   void fitsReportError(FILE * stream, int status) {
      if (status != 0) {
         fits_report_error(stream, status);
         throw std::runtime_error("ModelMap::trimExtensions(): "
                                  "cfitsio error.");
      }
   }

   void fitsResizeImage(const std::string & filename,
                        int bitpix, int naxis, 
                        std::vector<long> naxes) {
      fitsfile * fptr;
      int status(0);
      fits_open_file(&fptr, filename.c_str(), READWRITE, &status);
      fitsReportError(stderr, status);

      fits_resize_img(fptr, bitpix, naxis, &naxes[0], &status);
      fitsReportError(stderr, status);
      
      fits_close_file(fptr, &status);
      fitsReportError(stderr, status);
   }

   void toUpper(std::string & name) {
      for (std::string::iterator it = name.begin(); it != name.end(); ++it) {
         *it = std::toupper(*it);
      }
   }
} // anonymous namespace

namespace Likelihood {

ModelMap::ModelMap(BinnedLikelihood & logLike,
                   const std::vector<float> * model_map) 
   : m_logLike(logLike) {
   if (model_map == 0) {
      m_logLike.computeModelMap(m_outmap);
   } else {
      m_outmap = *model_map;
   }
}

void ModelMap::writeOutputMap(const std::string & outfile,
			      std::string outtype){
  ::toUpper(outtype);

  switch ( m_logLike.countsMap().projection().method() ){
  case astro::ProjBase::WCS:
    writeOutputMap_wcs(outfile,static_cast<const CountsMap&>(m_logLike.countsMap()),outtype);
    return;
  case astro::ProjBase::HEALPIX:
    writeOutputMap_healpix(outfile,static_cast<const CountsMapHealpix&>(m_logLike.countsMap()),outtype);
    return;
  default:
    break;
  } 
  std::ostringstream message;
  message << "Could not recognize projection method for CountsMap." ;
  throw std::runtime_error(message.str());
  return;
}

void ModelMap::writeOutputMap_healpix(const std::string & outfile,
				      const CountsMapHealpix& cmap,
				      std::string outtype) {
  
  ::toUpper(outtype);
   if (outtype != "CMAP" && outtype != "CCUBE") {
     std::ostringstream message;
      message << "Invalid output file type for model map, '" 
              << outtype << "'.\n"
              << "Only 'CMAP' and 'CCUBE' are allowed.";
      throw std::runtime_error(message.str());
   }
   std::string infile(cmap.filename());
   if (outtype == "CCUBE") {
     CountsMapHealpix outmap(cmap);
     outmap.setImage(m_outmap);
     outmap.writeOutput("gtmodel",outfile);
     return;
   } else {
     // use the c'tor to automatically do the sum over energies
     int nEBins = std::max( cmap.energies().size() -1, size_t(1));
     CountsMapHealpix outmap(cmap,0,nEBins);
     outmap.setImage(m_outmap);
     outmap.writeOutput("gtmodel",outfile);
   } 
   // Add the GTIs to the file
   dataSubselector::Cuts my_cuts(infile, "", false);
   my_cuts.writeGtiExtension(outfile);
   return;
}



void ModelMap::writeOutputMap_wcs(const std::string & outfile,
				      const CountsMap& cmap,
				      std::string outtype) {
 
   std::clock_t start = std::clock();
   
   ::toUpper(outtype);
   if (outtype != "CMAP" && outtype != "CCUBE") {
      std::ostringstream message;
      message << "Invalid output file type for model map, '" 
              << outtype << "'.\n"
              << "Only 'CMAP' and 'CCUBE' are allowed.";
      throw std::runtime_error(message.str());
   }

   if (outtype == "CMAP") {
      // Sum up the image planes over the different energy bands.
      size_t image_size(cmap.imageDimension(0)*
                        cmap.imageDimension(1));
      size_t nee(cmap.imageDimension(2));
      for (size_t k(1); k < nee; k++) {
         for (size_t j(0); j < image_size; j++) {
            size_t indx(k*image_size + j);
            m_outmap[j] += m_outmap[indx];
         }
      }
   }

   cmap.writeOutput("gtmodel", outfile);

   std::string infile(cmap.filename());
   const tip::Image * orig_image = tip::IFileSvc::instance().readImage(infile, "");
   const tip::Header & orig_header = orig_image->getHeader();

   std::vector<long> new_dims;
   typedef std::vector<tip::PixOrd_t> DimCont_t;
   new_dims.push_back(cmap.naxis1());
   new_dims.push_back(cmap.naxis2());
   if (outtype == "CCUBE") {
      new_dims.push_back(cmap.num_ebins());
   }

   ::fitsResizeImage(outfile, -32, new_dims.size(), new_dims);

   std::auto_ptr<tip::Image> output_image(tip::IFileSvc::instance().editImage(outfile, ""));
   tip::Header & output_header = output_image->getHeader();
   output_image->set(m_outmap);
   
   for ( tip::Header::ConstIterator itr = orig_header.begin(); itr != orig_header.end(); itr++ ) {
     // Don't copy the NAXIS* keywords
     if ( itr->getName().find("NAXIS") != std::string::npos ) {
       continue;
     }
     output_header.append(*itr);
   }
   if (outtype == "CMAP") {
     output_header.erase("CTYPE3");
     output_header.erase("CRPIX3");
     output_header.erase("CRVAL3");
     output_header.erase("CDELT3");
     output_header.erase("CUNIT3");
   }

}


} // namespace Likelihood
