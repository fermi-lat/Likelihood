/**
 * @file SourceMap.cxx
 * @brief Spatial distribution of a source folded through the instrument
 *        response.
 * @author J. Chiang
 *
 * $Header$
 */

#include "Likelihood/CountsMap.h"
#include "Likelihood/ExposureCube.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/Source.h"
#include "Likelihood/SourceMap.h"

namespace Likelihood {

SourceMap::SourceMap(Source * src, const CountsMap & dataMap) 
   : m_name(src->getName()) {
   std::vector<Pixel> pixels;
   dataMap.getPixels(pixels);
   
   std::vector<double> energies;
   dataMap.getAxisVector(2, energies);

// @todo Generalize this for possible future event types, perhaps
// making EventTypes a Singleton class in irfsInterface. For now, 0 =
// Front, 1 = Back.

   for (int evtType = 0; evtType < 2; evtType++) {
      std::vector<double> model;
      std::vector<double>::const_iterator energy = energies.begin();
      for ( ; energy != energies.end(); ++energy) {
         std::vector<Pixel>::const_iterator pixel = pixels.begin();
         for ( ; pixel != pixels.end(); ++pixel) {
            Aeff aeff(src, pixel->dir(), *energy, evtType);
            model.push_back(ExposureCube::instance()->value(pixel->dir(), 
                                                            aeff));
         }
      }
      m_models.push_back(model);
   }
}

double SourceMap::Aeff::operator()(double costheta) const {
   double inclination = acos(costheta);
   static double phi(0);
   return ResponseFunctions::totalResponse(inclination, phi, m_energy,
                                           m_energy, m_separation, m_type);
}

} // namespace Likelihood
