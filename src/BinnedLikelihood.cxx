/**
 * @file BinnedLikelihood.cxx
 * @brief First cut at a binned likelihood implementation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/BinnedLikelihood.cxx,v 1.8 2004/09/22 20:05:37 jchiang Exp $
 */

#include <stdexcept>

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Image.h"

#include "Likelihood/BinnedLikelihood.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/ExposureCube.h"
#include "Likelihood/SourceMap.h"
#include "Likelihood/SourceModel.h"

namespace Likelihood {

BinnedLikelihood::BinnedLikelihood(const CountsMap & dataMap,
                                   const std::string & srcMapsFile)
   : m_dataMap(dataMap), m_modelIsCurrent(false), m_srcMapsFile(srcMapsFile) {
   dataMap.getPixels(m_pixels);
   dataMap.getAxisVector(2, m_energies);
}

double BinnedLikelihood::value(optimizers::Arg &dummy) const {

   (void)(dummy);

   computeModelMap(m_energies, m_model);
   m_modelIsCurrent = true;
   
   const std::vector<double> & data = m_dataMap.data();
   double my_value(0);

   for (unsigned int j = 0; j < data.size(); j++) {
      if (m_model[j] > 0) {
         my_value += data[j]*log(m_model[j]) - m_model[j];
      }
   }
   
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
   if (!m_modelIsCurrent) {
      computeModelMap(m_energies, m_model);
      m_modelIsCurrent = true;
   }
   const std::vector<double> & data = m_dataMap.data();
   for (unsigned int j = 0; j < m_pixels.size(); j++) {
      for (unsigned int k = 0; k < m_energies.size()-1; k++) {
         unsigned long imin = k*m_pixels.size() + j;
         unsigned long imax = (k+1)*m_pixels.size() + j;
         if (m_model[imin] > 0) {
            long iparam(0);
            const std::map<std::string, Source *> & my_srcs = sources();
            std::map<std::string, Source *>::const_iterator src;
            for (src = my_srcs.begin(); src != my_srcs.end(); ++src) {
               std::string srcName = src->second->getName();
               const std::vector<double> & model = sourceMap(srcName)->model();
// This implementation assumes that only the spectral part of each
// Source has free parameters.
               Source::FuncMap srcFuncs = src->second->getSrcFuncs();
               Source::FuncMap::const_iterator func = srcFuncs.begin();
               for ( ; func != srcFuncs.end(); ++func) {
                  std::vector<std::string> paramNames;
                  func->second->getFreeParamNames(paramNames);
                  for (unsigned int i = 0; i < paramNames.size(); i++) {
                     double my_deriv = m_pixels[j].solidAngle()
                        *src->second->pixelCountsDeriv(m_energies[k], 
                                                       m_energies[k+1],
                                                       model[imin], 
                                                       model[imax],
                                                       paramNames[i]);
                     derivs.at(iparam) += (data[imin]/m_model[imin] - 1.)
                        *my_deriv;
                     iparam++;
                  }
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

   std::vector<double> energies;
   m_dataMap.getAxisVector(2, energies);

   std::vector<double> map;
   computeModelMap(energies, map);

   CountsMap * modelMap = new CountsMap(m_dataMap);
// @bug This will not work properly until
// tip::FitsExtensionManager::setImageDimensions is fixed.
   modelMap->setImage(map);
   return modelMap;
}

void BinnedLikelihood::computeModelMap(const std::vector<double> &energies,
                                       std::vector<double> & modelMap) const {
   modelMap.clear();
   modelMap.resize(m_pixels.size()*(energies.size()-1), 0);
   for (unsigned int k = 0; k < energies.size()-1; k++) {
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
            modelMap[imin] += m_pixels[j].solidAngle()
               *src->pixelCounts(energies[k], energies[k+1],
                                 model[imin], model[imax]);
         }
      }
   }
}

void BinnedLikelihood::createSourceMaps() {
   if (m_srcMapsFile != "") {
// read in the SourceMaps from the file
   } else {
      std::vector<std::string> srcNames;
      getSrcNames(srcNames);
      std::vector<std::string>::const_iterator name = srcNames.begin();
      for ( ; name != srcNames.end(); ++name) {
         Source * src = getSource(*name);
         m_srcMaps[*name] = new SourceMap(src, m_dataMap);
      }
   }
}

void BinnedLikelihood::saveSourceMaps(std::string filename) const {

   char * root_dir = std::getenv("LIKELIHOODROOT");
   if (root_dir == 0) {
      throw std::runtime_error("LIKELIHOODROOT not set.");
   }

   std::string templateFile = std::string(root_dir) 
      + "/data/LatCountsMapTemplate";

   if (filename == "") {
      if (m_srcMapsFile == "") {
         throw std::runtime_error("BinnedLikelihood::saveSourceMaps: "
                                  + std::string("filename not specified."));
      }
      filename = m_srcMapsFile;
   }

   tip::IFileSvc::instance().createFile(filename, templateFile);

   tip::Image * image = tip::IFileSvc::instance().editImage(filename, "");

   long image_dims[] = {m_dataMap.imageDimension(0), 
                        m_dataMap.imageDimension(1),
                        m_energies.size()};
   setImageDimensions(image, image_dims);

   tip::Header & header = image->getHeader();
   header["NAXIS1"].set(image_dims[0]);
   header["NAXIS2"].set(image_dims[1]);
   header["NAXIS3"].set(image_dims[2]);

   m_dataMap.setKeywords(header);

   std::vector<std::string> srcNames;
   getSrcNames(srcNames);
   const std::vector<double> & model 
      = m_srcMaps.find(srcNames[0])->second->model();
   std::vector<float> image_data(model.size());
   std::copy(model.begin(), model.end(), image_data.begin());
   image->set(image_data);

   delete image;
}

void BinnedLikelihood::setImageDimensions(tip::Image * image, 
                                          long * dimensions) const {
   typedef std::vector<tip::PixOrd_t> DimCont_t;
   DimCont_t dims = image->getImageDimensions();
   DimCont_t::size_type num_dims = dims.size();
   if (3 != num_dims) {
      throw std::runtime_error("BinnedLikelihood::setImageDimension: "
                               + std::string("cannot write a SourceMap file ")
                               + "to an image which is not 3D");
   }
   for (DimCont_t::size_type index = 0; index != num_dims; ++index) {
      dims[index] = dimensions[index];
   }
   image->setImageDimensions(dims);
}


} // namespace Likelihood
