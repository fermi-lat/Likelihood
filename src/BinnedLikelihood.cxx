/**
 * @file BinnedLikelihood.cxx
 * @brief Photon events are binned in sky direction and energy.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/BinnedLikelihood.cxx,v 1.12 2004/09/25 16:38:08 jchiang Exp $
 */

#include <memory>
#include <stdexcept>

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Image.h"

#include "st_facilities/Util.h"

#include "Likelihood/BinnedLikelihood.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/ExposureCube.h"
#include "Likelihood/SourceMap.h"
#include "Likelihood/SourceModel.h"

namespace Likelihood {

#include "fitsio.h"

BinnedLikelihood::BinnedLikelihood(const CountsMap & dataMap,
                                   const std::string & srcMapsFile)
   : m_dataMap(dataMap), m_modelIsCurrent(false), m_srcMapsFile(srcMapsFile) {
   dataMap.getPixels(m_pixels);
   dataMap.getAxisVector(2, m_energies);
   identifyFilledPixels();
}

double BinnedLikelihood::value(optimizers::Arg &dummy) const {
   (void)(dummy);

   double npred;
   computeModelMap(npred);

   const std::vector<double> & data = m_dataMap.data();
   double my_value(0);

   for (unsigned int i = 0; i < m_filledPixels.size(); i++) {
      if (m_model[i] > 0) {
         unsigned int j = m_filledPixels[i];
         my_value += data[j]*log(m_model[i]);
      }
   }
   my_value -= npred;

   return my_value;
}

std::vector<double>::const_iterator 
BinnedLikelihood::setParamValues_(std::vector<double>::const_iterator it) {
   m_modelIsCurrent = false;
   return SourceModel::setParamValues_(it);
}

std::vector<double>::const_iterator 
BinnedLikelihood::setFreeParamValues_(std::vector<double>::const_iterator it) {
   m_modelIsCurrent = false;
   return SourceModel::setFreeParamValues_(it);
}

void BinnedLikelihood::getFreeDerivs(std::vector<double> & derivs) const {
   int nparams(getNumFreeParams());
   derivs.resize(nparams, 0);
   double npred;
   if (!m_modelIsCurrent) {
      computeModelMap(npred);
   }
   const std::vector<double> & data = m_dataMap.data();
   for (unsigned int j = 0; j < m_filledPixels.size(); j++) {
      unsigned long imin = m_filledPixels[j];
      unsigned long imax = imin + m_pixels.size();
      unsigned long k = imin/m_pixels.size();
      if (m_model[j] > 0) {
         long iparam(0);
         const std::map<std::string, Source *> & my_srcs = sources();
         std::map<std::string, Source *>::const_iterator src;
         for (src = my_srcs.begin(); src != my_srcs.end(); ++src) {
            std::string srcName = src->second->getName();
            const SourceMap * srcMap = sourceMap(srcName);
            const std::vector<double> & model = srcMap->model();
// This implementation assumes that only the spectral part of each
// Source has free parameters.
            Source::FuncMap srcFuncs = src->second->getSrcFuncs();
            Source::FuncMap::const_iterator func = srcFuncs.begin();
            for ( ; func != srcFuncs.end(); ++func) {
               std::vector<std::string> paramNames;
               func->second->getFreeParamNames(paramNames);
               for (unsigned int i = 0; i < paramNames.size(); i++) {
                  double my_deriv = 
                     src->second->pixelCountsDeriv(m_energies[k], 
                                                   m_energies[k+1],
                                                   model[imin], 
                                                   model[imax],
                                                   paramNames[i]);
                  derivs.at(iparam) += data[imin]/m_model[j]*my_deriv;
                  if (j == 0) {
                     const std::vector<double> & npreds = srcMap->npreds();
                     for (unsigned int k = 0; k < m_energies.size()-1; k++) {
                        derivs.at(iparam) -= 
                           src->second->pixelCountsDeriv(m_energies[k], 
                                                         m_energies[k+1],
                                                         npreds[k], 
                                                         npreds[k+1],
                                                         paramNames[i]);
                     }
                  }
                  iparam++;
               }
            }
         }
      }
   }
}

CountsMap * BinnedLikelihood::createCountsMap() const {
   if (ExposureCube::instance() == 0) {
      std::runtime_error("BinnedLikelihood::createCountsMap:\n"
                         + std::string("Exposure cube not available."));
   }
   std::vector<double> map;
   computeModelMap(map);

   CountsMap * modelMap = new CountsMap(m_dataMap);
         
// @bug This will not work properly until
// tip::FitsExtensionManager::setImageDimensions is fixed.
   modelMap->setImage(map);
   return modelMap;
}

void BinnedLikelihood::computeModelMap(std::vector<double> & modelMap) const {
   modelMap.clear();
   modelMap.resize(m_pixels.size()*(m_energies.size()-1), 0);
   for (unsigned int k = 0; k < m_energies.size()-1; k++) {
      for (unsigned int j = 0; j < m_pixels.size(); j++) {
         unsigned long imin = k*m_pixels.size() + j;
         unsigned long imax = (k+1)*m_pixels.size() + j;
         std::map<std::string, SourceMap *>::const_iterator srcIt
            = m_srcMaps.begin();
         for ( ; srcIt != m_srcMaps.end(); ++srcIt) {
            const std::string & name = srcIt->first;
            const SourceMap * srcMap = srcIt->second;
            const Source * src 
               = const_cast<BinnedLikelihood *>(this)->getSource(name);
            const std::vector<double> & model = srcMap->model();
            modelMap[imin] += src->pixelCounts(m_energies[k], m_energies[k+1],
                                               model[imin], model[imax]);
         }
      }
   }
}

void BinnedLikelihood::computeModelMap(double & npred) const {
   m_model.clear();
   m_model.resize(m_filledPixels.size(), 0);
   npred = 0;
   std::map<std::string, SourceMap *>::const_iterator srcIt
      = m_srcMaps.begin();
   for ( ; srcIt != m_srcMaps.end(); ++srcIt) {
      const std::string & name = srcIt->first;
      const SourceMap * srcMap = srcIt->second;
      const Source * src 
         = const_cast<BinnedLikelihood *>(this)->getSource(name);
      if (src) {
         const std::vector<double> & model = srcMap->model();
         for (unsigned int i = 0; i < m_filledPixels.size(); i++) {
            unsigned long imin = m_filledPixels[i];
            unsigned long imax = imin + m_pixels.size();
            unsigned long k = imin/m_pixels.size();
            m_model.at(i) += src->pixelCounts(m_energies[k], m_energies[k+1],
                                              model[imin], model[imax]);
         }
         const std::vector<double> & npreds = srcMap->npreds();
         for (unsigned int k = 0; k < m_energies.size()-1; k++) {
            npred += src->pixelCounts(m_energies[k], m_energies[k+1],
                                      npreds[k], npreds[k+1]);
         }
      }
   }
   m_modelIsCurrent = true;
}

void BinnedLikelihood::createSourceMaps() {
   std::vector<std::string> srcNames;
   getSrcNames(srcNames);
   std::vector<std::string>::const_iterator name = srcNames.begin();
   for ( ; name != srcNames.end(); ++name) {
      Source * src = getSource(*name);
      m_srcMaps[*name] = new SourceMap(src, m_dataMap);
   }
}

void BinnedLikelihood::readSourceMaps(std::string filename) {
   if (filename == "") {
      if (m_srcMapsFile == "") {
         throw std::runtime_error("BinnedLikelihood::readSourceMaps: " +
                                  std::string("need to specify a SourceMaps ") 
                                  + "file.");
      }
      filename = m_srcMapsFile;
   }
   std::vector<std::string> srcNames;
   getSrcNames(srcNames);
   std::vector<std::string>::const_iterator name = srcNames.begin();
   for ( ; name != srcNames.end(); ++name) {
      if (sourceMapExists(*name, m_srcMapsFile)) {
// Replace an existing map in memory.
         if (m_srcMaps.count(*name)) {
            delete m_srcMaps[*name];
         }
/// @todo Confirm that the CountsMap in m_srcMapsFile matches m_dataMap.
         m_srcMaps[*name] = new SourceMap(m_srcMapsFile, *name);
      }
   }
// Generate needed maps that not available in the file.
   for (name = srcNames.begin() ; name != srcNames.end(); ++name) {
      if (!m_srcMaps.count(*name)) {
         Source * src = getSource(*name);
         m_srcMaps[*name] = new SourceMap(src, m_dataMap);
      }
   }
}

void BinnedLikelihood::saveSourceMaps(const std::string & filename) {
   if (filename != "") {
      m_srcMapsFile = filename;
   }
   if (!st_facilities::Util::fileExists(filename)) {
      m_dataMap.writeOutput("BinnedLikelihood", m_srcMapsFile);
   }
   std::vector<std::string> srcNames;
   getSrcNames(srcNames);
   for (unsigned int i = 0; i < srcNames.size(); i++) {
      if (sourceMapExists(srcNames[i], m_srcMapsFile)) {
         replaceSourceMap(srcNames[i], m_srcMapsFile);
      } else {
         addSourceMap(srcNames[i], m_srcMapsFile);
      }
   }
}

bool BinnedLikelihood::sourceMapExists(const std::string & srcName,
                                       const std::string & fitsFile) const {
   try {
      std::auto_ptr<const tip::Image> 
         image(tip::IFileSvc::instance().readImage(fitsFile, srcName));
   } catch (tip::TipException & eObj) {
      return false;
   }
   return true;
}

void BinnedLikelihood::replaceSourceMap(const std::string & srcName,
                                        const std::string & fitsFile) const {
   (void)(srcName);
   (void)(fitsFile);
/// @todo Implement replaceSourceMap method.
}

void BinnedLikelihood::addSourceMap(const std::string & srcName,
                                    const std::string & fitsFile) const {
   if (!m_srcMaps.count(srcName)) {
      throw std::runtime_error("BinnedLikelihood::addSourceMap: " +
                               std::string("Source ") + srcName 
                               + " is not available.");
   }

   fitsfile * fptr;
   int status(0);

   fits_open_file(&fptr, fitsFile.c_str(), READWRITE, &status);
   fitsReportError(stderr, status);

   long naxes[] = {m_dataMap.imageDimension(0),
                   m_dataMap.imageDimension(1),
                   m_energies.size()};
   fits_create_img(fptr, DOUBLE_IMG, 3, naxes, &status);
   fitsReportError(stderr, status);
   
   long group(0);
   const std::vector<double> & data = m_srcMaps.find(srcName)->second->model();
   fits_write_3d_dbl(fptr, group, naxes[0], naxes[1], naxes[0], naxes[1],
                     naxes[2], const_cast<double *>(&data[0]), &status);
   fitsReportError(stderr, status);

   fits_update_key(fptr, TSTRING, "EXTNAME", 
                   const_cast<char *>(srcName.c_str()), "SourceMap name",
                   &status);
   fitsReportError(stderr, status);

   fits_close_file(fptr, &status);
   fitsReportError(stderr, status);
}

void BinnedLikelihood::fitsReportError(FILE *stream, int status) const {
   if (status != 0) {
      fits_report_error(stream, status);
      throw std::runtime_error("BinnedLikelihood::addSourceMap: "
                               + std::string("cfitsio error."));
   }
}

void BinnedLikelihood::setImageDimensions(tip::Image * image, 
                                          long * dimensions) const {
   typedef std::vector<tip::PixOrd_t> DimCont_t;
   DimCont_t dims = image->getImageDimensions();
   DimCont_t::size_type num_dims = dims.size();
   if (3 != num_dims) {
      throw std::runtime_error("BinnedLikelihood::setImageDimensions: "
                               + std::string("Cannot write a SourceMap file ")
                               + "to an image which is not 3D");
   }
   for (DimCont_t::size_type index = 0; index != num_dims; ++index) {
      dims[index] = dimensions[index];
   }
   image->setImageDimensions(dims);
}

void BinnedLikelihood::identifyFilledPixels() {
   const std::vector<double> & data = m_dataMap.data();
   m_filledPixels.clear();
   for (unsigned int i = 0; i < data.size(); i++) {
      if (data.at(i) > 0) {
         m_filledPixels.push_back(i);
      }
   }
}

} // namespace Likelihood
