/**
 * @file ModelMap.cxx
 * @brief Standalone class for compute a model map using a 
 * BinnedLikelihood object.
 */

#include <memory>
#include <sstream>
#include <stdexcept>

#include "fitsio.h"

#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Table.h"

#include "st_facilities/FitsUtil.h"

#include "dataSubselector/Cuts.h"

#include "Likelihood/BinnedLikelihood.h"

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
                              std::string outtype) {
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
      size_t image_size(m_logLike.countsMap().imageDimension(0)*
                        m_logLike.countsMap().imageDimension(1));
      size_t nee(m_logLike.countsMap().imageDimension(2));
      for (size_t k(1); k < nee; k++) {
         for (size_t j(0); j < image_size; j++) {
            size_t indx(k*image_size + j);
            m_outmap[j] += m_outmap[indx];
         }
      }
   }
   std::string infile(m_logLike.countsMap().filename());
   bool clobber;
   st_facilities::FitsUtil::fcopy(infile, outfile, "", "", clobber=true);
   
   std::vector<long> new_dims;
   const tip::Image * my_image = 
      tip::IFileSvc::instance().readImage(outfile, "");
   typedef std::vector<tip::PixOrd_t> DimCont_t;
   DimCont_t dims = my_image->getImageDimensions();
   new_dims.push_back(dims.at(0));
   new_dims.push_back(dims.at(1));
   if (outtype == "CCUBE") {
      new_dims.push_back(dims.at(2));
   }
   delete my_image;
   ::fitsResizeImage(outfile, -32, new_dims.size(), new_dims);

   std::auto_ptr<tip::Image>
      output_image(tip::IFileSvc::instance().editImage(outfile, ""));
   output_image->set(m_outmap);
   dataSubselector::Cuts my_cuts(infile, "", false);
   my_cuts.writeGtiExtension(outfile);
   trimExtensions(outfile, outtype);
}

void ModelMap::trimExtensions(const std::string & outfile,
                                 const std::string & outtype) {
   tip::FileSummary hdus;
   tip::IFileSvc::instance().getFileSummary(outfile, hdus);
   size_t nhdus = hdus.size();
   fitsfile * fptr;
   int status(0), hdutype;
   fits_open_file(&fptr, outfile.c_str(), READWRITE, &status);
   ::fitsReportError(stderr, status);
   for (size_t i = 1; i < nhdus; i++) {
      std::string hduname = hdus.at(i).getExtId();
      if (hduname != "GTI" && (hduname != "EBOUNDS" || outtype == "CMAP")) {
         fits_movnam_hdu(fptr, ANY_HDU, const_cast<char *>(hduname.c_str()), 
                         0, &status);
         ::fitsReportError(stderr, status);
         fits_delete_hdu(fptr, &hdutype, &status);
         ::fitsReportError(stderr, status);
      }
   }
   if (outtype == "CMAP") {
// Delete axis 3 keywords in PRIMARY HDU.
      fits_movabs_hdu(fptr, 1, &hdutype, &status);
      ::fitsReportError(stderr, status);

      fits_delete_key(fptr, "CTYPE3", &status);
      ::fitsReportError(stderr, status);
      fits_delete_key(fptr, "CRPIX3", &status);
      ::fitsReportError(stderr, status);
      fits_delete_key(fptr, "CRVAL3", &status);
      ::fitsReportError(stderr, status);
      fits_delete_key(fptr, "CDELT3", &status);
      ::fitsReportError(stderr, status);
      fits_delete_key(fptr, "CUNIT3", &status);
      ::fitsReportError(stderr, status);
   }

// Update creator keyword.
   char * keyname = "CREATOR";
   char * creator = "gtmodel";
   char * description = "Software creating file";
   fits_update_key(fptr, TSTRING, keyname, creator, description, &status);
   ::fitsReportError(stderr, status);

   fits_close_file(fptr, &status);
   ::fitsReportError(stderr, status);

   st_facilities::FitsUtil::writeChecksums(outfile);
}

} // namespace Likelihood
