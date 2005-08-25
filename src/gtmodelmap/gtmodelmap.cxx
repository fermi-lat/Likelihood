/**
 * @file gtmodelmap.cxx
 * @brief Compute a model counts map based on binned likelihood fits.
 * @author J. Chiang
 *
 * $Header$
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

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "st_facilities/FitsUtil.h"

#include "optimizers/dArg.h"
#include "optimizers/FunctionFactory.h"

#include "Likelihood/BandFunction.h"
#include "Likelihood/LogParabola.h"
#include "Likelihood/MapCubeFunction.h"
#include "Likelihood/PowerLaw2.h"
#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/SpatialMap.h"

#include "fitsio.h"

XERCES_CPP_NAMESPACE_USE

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
         throw std::runtime_error("gtmodelmap::trimExtensions(): "
                                  "cfitsio error.");
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
                m_pars(st_app::StApp::getParGroup("gtmodelmap")),
                m_funcFactory(0) {}
   virtual ~ModelMap() throw() {
      try {
         delete m_funcFactory;
      } catch (std::exception & eObj) {
         std::cout << eObj.what() << std::endl;
      } catch (...) {
      }
   }
   virtual void run();
   virtual void banner() const {};

private:

   st_app::AppParGroup & m_pars;

   optimizers::FunctionFactory * m_funcFactory;

   std::map<std::string, optimizers::Function *> m_spectra;

   std::vector<double> m_emins;
   std::vector<double> m_emaxs;

   std::vector<float> m_srcmap;
   std::vector<float> m_outmap;

   void prepareFunctionFactory();
   void readSpectra();
   void readEnergyBounds();
   void sumOutputMap();
   void writeOutputMap();
   void trimExtensions();

   void getMap(const std::string & srcName);
                
};

st_app::StAppFactory<ModelMap> myAppFactory("gtmodelmap");

void ModelMap::run() {
   m_pars.Prompt();
   m_pars.Save();
   prepareFunctionFactory();
   readSpectra();
   readEnergyBounds();
   sumOutputMap();
   writeOutputMap();
   trimExtensions();
}

void ModelMap::writeOutputMap() {
   std::string infile = m_pars["srcmaps"];
   std::string outfile = m_pars["outfile"];
   bool clobber = m_pars["clobber"];
   tip::IFileSvc::instance().createFile(outfile, infile, clobber);
   std::auto_ptr<tip::Image>
      output_image(tip::IFileSvc::instance().editImage(outfile, ""));
   typedef std::vector<tip::PixOrd_t> DimCont_t;
   DimCont_t dims = output_image->getImageDimensions();
   dims.resize(2);
   output_image->setImageDimensions(dims);
   output_image->set(m_outmap);
}

void ModelMap::trimExtensions() {
   std::string outfile = m_pars["outfile"];
   tip::FileSummary hdus;
   tip::IFileSvc::instance().getFileSummary(outfile, hdus);
   size_t nhdus = hdus.size();
   fitsfile * fptr;
   int status(0), hdutype;
   fits_open_file(&fptr, outfile.c_str(), READWRITE, &status);
   ::fitsReportError(stderr, status);
   for (size_t i = 1; i < nhdus; i++) {
      std::string hduname = hdus.at(i).getExtId();
      if (hduname != "GTI") {
         fits_movnam_hdu(fptr, ANY_HDU, const_cast<char *>(hduname.c_str()), 
                         0, &status);
         ::fitsReportError(stderr, status);
         fits_delete_hdu(fptr, &hdutype, &status);
         ::fitsReportError(stderr, status);
      }
   }
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

// update creator keyword
   char * creator = "gtmodelmap";
   fits_update_key(fptr, TSTRING, "CREATOR", creator,
                   "Software creating file", &status);
   ::fitsReportError(stderr, status);

   fits_close_file(fptr, &status);
   ::fitsReportError(stderr, status);

   st_facilities::FitsUtil::writeChecksums(outfile);
}

void ModelMap::sumOutputMap() {
   std::map<std::string, optimizers::Function *>::iterator it;
   for (it = m_spectra.begin(); it != m_spectra.end(); ++it) {
      std::string srcName = it->first;
      try {
         getMap(srcName);
      } catch (tip::TipException & eObj) {
         std::cout << "Cannot read source map for model component "
                   << srcName << ". Skipping it." << std::endl;
      }
      if (it == m_spectra.begin()) {
         m_outmap.resize(m_srcmap.size(), 0);
      }
      ::Spectrum spec(it->second);
      size_t image_size = m_outmap.size()/(m_emins.size() + 1);
      for (unsigned int k = 0; k < m_emins.size(); k++) {
         double & emin = m_emins.at(k);
         double & emax = m_emaxs.at(k);
         double de = (emax - emin)/2.;
         for (unsigned int i = 0; i < image_size; i++) {
            size_t j0 = k*image_size + i;
            size_t j1 = j0 + image_size;
            m_outmap.at(i) += (m_srcmap.at(j0)*spec(emin) +
                               m_srcmap.at(j1)*spec(emax))*de;
         }
      }
   }
}

void ModelMap::getMap(const std::string & srcName) {
   std::string srcMaps_file = m_pars["srcmaps"];
   std::auto_ptr<const tip::Image> 
      image(tip::IFileSvc::instance().readImage(srcMaps_file, srcName));
   image->get(m_srcmap);
}

void ModelMap::readSpectra() {
   std::string model_file = m_pars["source_model_file"];

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
   bool makeClone(false);
   m_funcFactory = new optimizers::FunctionFactory();
   m_funcFactory->addFunc("SkyDirFunction", new Likelihood::SkyDirFunction(),
                          makeClone);
   m_funcFactory->addFunc("SpatialMap", new Likelihood::SpatialMap(),
                          makeClone);
   m_funcFactory->addFunc("BandFunction", new Likelihood::BandFunction(),
                          makeClone);
   m_funcFactory->addFunc("LogParabola", new Likelihood::LogParabola(),
                          makeClone);
   m_funcFactory->addFunc("MapCubeFunction", new Likelihood::MapCubeFunction(),
                          makeClone);
   m_funcFactory->addFunc("PowerLaw2", new Likelihood::PowerLaw2(), makeClone);
}
