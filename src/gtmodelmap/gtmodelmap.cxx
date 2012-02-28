/**
 * @file gtmodelmap.cxx
 * @brief Compute a model counts map based on binned likelihood fits.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/gtmodelmap/gtmodelmap.cxx,v 1.33 2012/01/17 22:05:25 jchiang Exp $
 */

#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <xercesc/util/XercesDefs.hpp>

#include "xmlBase/Dom.h"
#include "xmlBase/XmlParser.h"

#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Table.h"

#include "st_stream/StreamFormatter.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "st_facilities/FitsUtil.h"

#include "dataSubselector/Cuts.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/BinnedLikelihood.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/SourceMap.h"

#include "fitsio.h"

namespace {
   void fitsReportError(FILE * stream, int status) {
      if (status != 0) {
         fits_report_error(stream, status);
         throw std::runtime_error("gtmodel::trimExtensions(): "
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
}

/**
 * @class ModelMap
 *
 * @brief Derived class of st_app::StApp for summing up source maps
 * with the spectal fit parameters from a binned likelihood analysis
 * applied.
 *
 * @author J. Chiang
 *
 */

class ModelMap : public st_app::StApp {

public:

   ModelMap() : st_app::StApp(),
                m_pars(st_app::StApp::getParGroup("gtmodel")),
                m_helper(0), m_logLike(0) {
      setVersion(s_cvs_id);
   }
   virtual ~ModelMap() throw() {
      try {
         delete m_logLike;
         delete m_dataMap;
         delete m_helper;
      } catch (std::exception & eObj) {
         std::cerr << eObj.what() << std::endl;
      } catch (...) {
      }
   }
   virtual void run();
   virtual void banner() const;

private:

   st_app::AppParGroup & m_pars;

   Likelihood::AppHelpers * m_helper;
   Likelihood::CountsMap * m_dataMap;
   Likelihood::BinnedLikelihood * m_logLike;

   std::vector<float> m_outmap;

   void computeModelMap();
   void writeOutputMap();
   void trimExtensions();

   static std::string s_cvs_id;
};

st_app::StAppFactory<ModelMap> myAppFactory("gtmodel");

std::string ModelMap::s_cvs_id("$Name:  $");

void ModelMap::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

void ModelMap::run() {
   m_pars.Prompt();
   m_pars.Save();
   computeModelMap();
   writeOutputMap();
   trimExtensions();
}

void ModelMap::computeModelMap() {
   m_helper = new Likelihood::AppHelpers(&m_pars, "BINNED");
   m_helper->observation().expCube().readExposureCube(m_pars["expcube"]);
   m_helper->setRoi(m_pars["srcmaps"], "", false);
   std::string cmapfile = m_pars["srcmaps"];
   m_dataMap = new Likelihood::CountsMap(cmapfile);
   bool computePointSources, apply_psf_corrections;
   bool performConvolution = m_pars["convol"];
   bool resample = m_pars["resample"];
   int resamp_factor = m_pars["rfactor"];
   double rfactor = static_cast<double>(resamp_factor);
   m_logLike = new Likelihood::BinnedLikelihood(*m_dataMap,
                                                m_helper->observation(),
                                                cmapfile, 
                                                computePointSources=true, 
                                                apply_psf_corrections=true,
                                                performConvolution,
                                                resample, rfactor);
   std::string bexpmap = m_pars["bexpmap"];
   Likelihood::AppHelpers::checkExposureMap(cmapfile, bexpmap);
   if (bexpmap != "none" && bexpmap != "") {
      Likelihood::SourceMap::setBinnedExposure(bexpmap);
   }
   bool requireExposure, addPointSources, loadMaps, createAllMaps;
   m_logLike->readXml(m_pars["srcmdl"], m_helper->funcFactory(),
                      requireExposure=false, addPointSources=true,
                      loadMaps=false, createAllMaps=true);
   m_logLike->computeModelMap(m_outmap);

   std::string outtype = m_pars["outtype"];
   if (outtype == "CMAP") {
      // Sum up the image planes over the different energy bands.
      size_t image_size(m_dataMap->imageDimension(0)*
                        m_dataMap->imageDimension(1));
      size_t nee(m_dataMap->imageDimension(2));
      for (size_t k(1); k < nee; k++) {
         for (size_t j(0); j < image_size; j++) {
            size_t indx(k*image_size + j);
            m_outmap[j] += m_outmap[indx];
         }
      }
   }
}

void ModelMap::writeOutputMap() {
   std::string infile = m_pars["srcmaps"];
   std::string outfile = m_pars["outfile"];
   bool clobber = m_pars["clobber"];
   std::string outtype = m_pars["outtype"];
   st_facilities::FitsUtil::fcopy(infile, outfile, "", "", clobber);
   
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
   fitsResizeImage(outfile, -32, new_dims.size(), new_dims);

   std::auto_ptr<tip::Image>
      output_image(tip::IFileSvc::instance().editImage(outfile, ""));
   output_image->set(m_outmap);
   dataSubselector::Cuts my_cuts(infile, "", false);
   my_cuts.writeGtiExtension(outfile);
}

void ModelMap::trimExtensions() {
   std::string outfile = m_pars["outfile"];
   std::string outtype = m_pars["outtype"];
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
   const char * keyname = "CREATOR";
   char * creator = "gtmodel";
   char * description = "Software creating file";
   fits_update_key(fptr, TSTRING, keyname, creator, description, &status);
   ::fitsReportError(stderr, status);

   fits_close_file(fptr, &status);
   ::fitsReportError(stderr, status);

   st_facilities::FitsUtil::writeChecksums(outfile);
}
