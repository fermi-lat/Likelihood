/**
 * @file SourceMap.cxx
 * @brief Spatial distribution of a source folded through the instrument
 *        response.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SourceMap.cxx,v 1.7 2004/09/28 04:32:25 jchiang Exp $
 */

#include <algorithm>
#include <memory>

#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Table.h"
#include "tip/tip_types.h"

#include "Likelihood/BinnedExposure.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/ExposureCube.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/Source.h"
#include "Likelihood/SourceMap.h"

namespace Likelihood {

MeanPsf * SourceMap::s_meanPsf(0);
BinnedExposure * SourceMap::s_binnedExposure(0);

SourceMap::SourceMap(Source * src, const CountsMap & dataMap) 
   : m_name(src->getName()), m_dataMap(dataMap) {
   std::vector<Pixel> pixels;
   dataMap.getPixels(pixels);
   
   std::vector<double> energies;
   dataMap.getAxisVector(2, energies);

   std::cerr << "Generating SourceMap for " << m_name;
   long npts = energies.size()*pixels.size();
   m_model.resize(npts, 0);
   long icount(0);

   m_npreds.resize(energies.size(), 0);

/// @todo Ensure the desired event types are correctly included in this
/// calculation.
   for (int evtType = 0; evtType < 2; evtType++) {
      std::vector<double>::const_iterator energy = energies.begin();
      for (int k = 0; energy != energies.end(); ++energy, k++) {
         std::vector<Pixel>::const_iterator pixel = pixels.begin();
         for (int j = 0; pixel != pixels.end(); ++pixel, j++) {
            unsigned long indx = k*pixels.size() + j;
            if ((icount % (npts/10)) == 0) std::cerr << ".";
            double value(0);
            if (dynamic_cast<PointSource *>(src) != 0) {
               Aeff aeff(src, pixel->dir(), *energy, evtType);
               value = ExposureCube::instance()->value(pixel->dir(), aeff)
                  *pixel->solidAngle();
            } else if ( dynamic_cast<DiffuseSource *>(src) != 0 ) {
               value = sourceRegionIntegral(src, *pixel, *energy, evtType);
            }
            m_model.at(indx) += value;
            m_npreds.at(k) += value;
            icount++;
         }
      }
   }
   std::cerr << "!" << std::endl;
}

SourceMap::SourceMap(const std::string & sourceMapsFile,
                     const std::string & srcName) 
   : m_name(srcName), m_dataMap(CountsMap(sourceMapsFile)) {
   std::auto_ptr<const tip::Image> 
      image(tip::IFileSvc::instance().readImage(sourceMapsFile, srcName));
   std::vector<float> image_data;
   image->get(image_data);
   m_model.resize(image_data.size());
   std::copy(image_data.begin(), image_data.end(), m_model.begin());

   std::vector<Pixel> pixels;
   m_dataMap.getPixels(pixels);
   std::vector<double> energies;
   m_dataMap.getAxisVector(2, energies);

   m_npreds.resize(energies.size(), 0);
   for (unsigned int k = 0; k < energies.size(); k++) {
      std::vector<Pixel>::const_iterator pixel = pixels.begin();
      for (int j = 0; pixel != pixels.end(); ++pixel, j++) {
         unsigned long indx = k*pixels.size() + j;
         m_npreds.at(k) += m_model.at(indx);
      }
   }
}

double SourceMap::Aeff::operator()(double costheta) const {
   double inclination = acos(costheta)*180./M_PI;
   static double phi(0);
   return ResponseFunctions::totalResponse(inclination, phi, m_energy,
                                           m_energy, m_separation, m_type);

}

double SourceMap::sourceRegionIntegral(Source * src, const Pixel & pixel,
                                       double energy, int evtType) const {
   DiffuseSource * diffuseSrc = dynamic_cast<DiffuseSource *>(src);
   std::vector<double> energies;
   m_dataMap.getAxisVector(2, energies);
   if (s_meanPsf == 0 || s_binnedExposure == 0) {
      delete s_meanPsf;
      delete s_binnedExposure;
      double ra, dec;
      RoiCuts::instance()->getRaDec(ra, dec);
      s_meanPsf = new MeanPsf(ra, dec, energies);
      s_binnedExposure = new BinnedExposure(energies);
   }
   MeanPsf & psf = *s_meanPsf;
   BinnedExposure & exposure = *s_binnedExposure;
   return 0;
}

} // namespace Likelihood
