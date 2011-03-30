/**
 * @file gtmodelmap.cxx
 * @brief Compute a model counts map based on binned likelihood fits.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/gtmodelmap/gtmodelmap.cxx,v 1.27 2011/02/02 01:21:49 jchiang Exp $
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

#include "optimizers/dArg.h"
#include "optimizers/FunctionFactory.h"

#include "dataSubselector/Cuts.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/FileFunction.h"
#include "Likelihood/DMFitFunction.h"
#include "Likelihood/DMFitFunction2.h"
#include "Likelihood/Source.h"

#include "SourceMapRegistry.h"

#include "fitsio.h"

//XERCES_CPP_NAMESPACE_USE
using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;
using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;

namespace {
   class Spectrum {
   public:
      Spectrum(optimizers::Function * func) : m_func(func) {}
      double operator()(double energy) {
         optimizers::dArg eArg(energy);
         return m_func->value(eArg);
      }
   private:
      optimizers::Function * m_func;
   };

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

   void tolower(std::string & name) {
      for (std::string::iterator it = name.begin(); it != name.end(); ++it) {
         *it = std::tolower(*it);
      }
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
                m_funcFactory(0), m_srcmap(0) {
      setVersion(s_cvs_id);
   }
   virtual ~ModelMap() throw() {
      try {
         delete m_funcFactory;
      } catch (std::exception & eObj) {
         std::cerr << eObj.what() << std::endl;
      } catch (...) {
      }
   }
   virtual void run();
   virtual void banner() const;

private:

   st_app::AppParGroup & m_pars;

   optimizers::FunctionFactory * m_funcFactory;

   std::map<std::string, optimizers::Function *> m_spectra;

   SourceMapRegistry * m_registry;

   std::vector<double> m_emins;
   std::vector<double> m_emaxs;

   std::vector<float> * m_srcmap;
   std::vector<float> m_outmap;

   void prepareFunctionFactory();
   void readSpectra();
   void readEnergyBounds();
   void createRegistry();
   void sumOutputMap();
   void writeOutputMap();
   void trimExtensions();

   void getMap(const std::string & srcName);

   double pixelCounts(double emin, double emax, double y1, double y2) const;

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
   prepareFunctionFactory();
   readSpectra();
   readEnergyBounds();
   createRegistry();
   sumOutputMap();
   writeOutputMap();
   trimExtensions();
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
// delete axis 3 keywords in PRIMARY HDU
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

// update creator keyword
   char * creator = "gtmodel";
   fits_update_key(fptr, TSTRING, "CREATOR", creator,
                   "Software creating file", &status);
   ::fitsReportError(stderr, status);

   fits_close_file(fptr, &status);
   ::fitsReportError(stderr, status);

   st_facilities::FitsUtil::writeChecksums(outfile);
}

void ModelMap::sumOutputMap() {
   std::string outtype = m_pars["outtype"];
   st_stream::StreamFormatter formatter("gtmodel", "sumOutputMap", 2);
   for (size_t k(0); k < m_emins.size(); k++) {
      double & emin(m_emins.at(k));
      double & emax(m_emaxs.at(k));
      std::vector<double> map0;
      std::vector<double> map1;
      size_t image_size(0);
      std::map<std::string, optimizers::Function *>::iterator it;
      for (it = m_spectra.begin(); it != m_spectra.end(); ++it) {
         std::string srcName = it->first;
         try {
            getMap(srcName);
            if (image_size == 0) {
               m_outmap.resize(m_srcmap->size(), 0);
               image_size = m_outmap.size()/(m_emins.size() + 1);
               map0.resize(image_size, 0);
               map1.resize(image_size, 0);
            }
         } catch (tip::TipException &) {
            formatter.info() << "Cannot read source map for model component "
                             << srcName << ". Skipping it." << std::endl;
         }
         for (size_t i(0); i < image_size; i++) {
            size_t j0 = k*image_size + i;
            optimizers::dArg eminArg(emin);
            map0.at(i) += it->second->operator()(eminArg)*m_srcmap->at(j0);

            size_t j1 = j0 + image_size;
            optimizers::dArg emaxArg(emax);
            map1.at(i) += it->second->operator()(emaxArg)*m_srcmap->at(j1);
         }
      }
      for (size_t i(0); i < image_size; i++) {
         if (outtype == "CMAP") {
            m_outmap.at(i) += pixelCounts(emin, emax, map0.at(i), map1.at(i));
         } else {
            size_t indx = i + k*image_size;
            m_outmap.at(indx) = pixelCounts(emin, emax, map0.at(i), map1.at(i));
         }
      }
   }
}

double ModelMap::pixelCounts(double emin, double emax,
                             double y1, double y2) const {
   if (::getenv("USE_LINEAR_QUADRATURE")) {
      return (y1 + y2)*(emax - emin)/2.;
   }
   return (y1*emin + y2*emax)/2.*std::log(emax/emin);
}

void ModelMap::getMap(const std::string & srcName) {
   if (m_registry) {
      m_srcmap = 
         const_cast<std::vector<float> *>(&m_registry->sourceMap(srcName));
   } else {
      std::string srcMaps_file = m_pars["srcmaps"];
      std::auto_ptr<const tip::Image> 
         image(tip::IFileSvc::instance().readImage(srcMaps_file, srcName));
      m_srcmap = new std::vector<float>(0);
      image->get(*m_srcmap);
   }
}

void ModelMap::readSpectra() {
   std::string model_file = m_pars["srcmdl"];

   xmlBase::XmlParser * parser = new xmlBase::XmlParser();
   DOMDocument * doc = parser->parse(model_file.c_str());

   if (doc == 0) {
      throw std::runtime_error("ModelMap::readSpectra(): "
                               "input model file not parsed successfully.");
   }

   DOMElement * source_model = doc->getDocumentElement();

   std::vector<DOMElement *> sources;
   xmlBase::Dom::getChildrenByTagName(source_model, "source", sources);

   for (unsigned int i=0; i < sources.size(); i++) {
      DOMElement * spectrum = 
         xmlBase::Dom::findFirstChildByName(sources.at(i), "spectrum");
      std::string type = xmlBase::Dom::getAttribute(spectrum, "type");
      optimizers::Function * func = m_funcFactory->create(type);
      func->setParams(spectrum);
      std::string name = xmlBase::Dom::getAttribute(sources.at(i), "name");

// If FileFunction or DMFitFunction[2], read in the data:
      if (type == "FileFunction") {
         std::string filename = xmlBase::Dom::getAttribute(spectrum, "file");
         dynamic_cast<Likelihood::FileFunction *>(func)->readFunction(filename);
      }
      if (type == "DMFitFunction") {
         std::string filename = xmlBase::Dom::getAttribute(spectrum, "file");
         dynamic_cast<Likelihood::DMFitFunction *>(func)->readFunction(filename);
      }
      if (type == "DMFitFunction2") {
         std::string filename = xmlBase::Dom::getAttribute(spectrum, "file");
         dynamic_cast<Likelihood::DMFitFunction2 *>(func)->readFunction(filename);
      }

      m_spectra[name] = func;
   }
}
                                                                 
void ModelMap::readEnergyBounds() {
   std::string srcMaps_file = m_pars["srcmaps"];
   const tip::Table * ebounds = 
      tip::IFileSvc::instance().readTable(srcMaps_file, "EBOUNDS");
   tip::Table::ConstIterator it = ebounds->begin();
   tip::ConstTableRecord & ebound = *it;
   m_emins.clear();
   m_emaxs.clear();
   double emin, emax;
   for ( ; it != ebounds->end(); ++it) {
      ebound["E_MIN"].get(emin);
      ebound["E_MAX"].get(emax);
      m_emins.push_back(emin/1e3);
      m_emaxs.push_back(emax/1e3);
   }
   delete ebounds;
}

void ModelMap::prepareFunctionFactory() {
   m_funcFactory = new optimizers::FunctionFactory();
   Likelihood::AppHelpers::addFunctionPrototypes(m_funcFactory);
}

void ModelMap::createRegistry() {
   std::string countsMap = m_pars["srcmaps"];
   std::string xmlFile = m_pars["srcmdl"];
   std::string irfs = m_pars["irfs"];
   std::string expCube = m_pars["expcube"];
   std::string binnedExpMap = m_pars["bexpmap"];
   bool performConvolution = m_pars["convol"];
   bool resample = m_pars["resample"];
   int factor = m_pars["rfactor"];
   double resamp_factor = static_cast<double>(factor);

   if (expCube != "none") {
      m_registry = new SourceMapRegistry(countsMap, xmlFile, irfs, expCube,
                                         binnedExpMap, *m_funcFactory,
                                         performConvolution, resample, 
                                         resamp_factor);
   } else {
      throw std::runtime_error("You must specify a livetime cube file.");
   }
}
