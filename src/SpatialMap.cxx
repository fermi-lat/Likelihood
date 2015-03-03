/** 
 * @file SpatialMap.cxx
 * @brief Implementation of Function object class that returns interpolated
 * image values of a FITS image file.
 * 
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SpatialMap.cxx,v 1.35 2013/09/04 05:30:45 jchiang Exp $
 *
 */

#include <algorithm>
#include <iostream>
#include <stdexcept>

#include "facilities/Util.h"

#include "st_stream/StreamFormatter.h"

#include "st_facilities/Util.h"

#include "Likelihood/Event.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/SkyDirArg.h"
#include "Likelihood/SpatialMap.h"

namespace Likelihood {

SpatialMap::SpatialMap() 
   : optimizers::Function("SpatialMap", 1, "Prefactor"), MapBase() {
   addParam("Prefactor", 1, false);
   setParamAlwaysFixed("Prefactor");
}

SpatialMap::SpatialMap(const std::string & fitsFile,
                       const std::string & extension) 
   : optimizers::Function("SpatialMap", 1, "Prefactor"), 
     MapBase(fitsFile, extension) {
   addParam("Prefactor", 1, false);
   setParamAlwaysFixed("Prefactor");
}

SpatialMap::SpatialMap(const SpatialMap & rhs) 
   : optimizers::Function(rhs), MapBase(rhs) {
}

SpatialMap & SpatialMap::operator=(const SpatialMap & rhs) {
   if (this != &rhs) {
      optimizers::Function::operator=(rhs);
      MapBase::operator=(rhs);
   }
   return *this;
}

SpatialMap::~SpatialMap() {
}

double SpatialMap::value(optimizers::Arg & arg) const {
   astro::SkyDir dir;
   dynamic_cast<SkyDirArg &>(arg).fetchValue(dir);
   return value(dir);
}

double SpatialMap::value(const astro::SkyDir & dir) const {
   double pref = m_parameter[0].getTrueValue();
   return pref*wcsmap().operator()(dir);
}

void SpatialMap::integrateSpatialDist(const std::vector<double> & energies,
                                      const ExposureMap & expmap,
                                      std::vector<double> & exposure) const {
   exposure.clear();
//    const std::vector< std::vector<double> > & 
//       solid_angles(wcsmap().solidAngles());
//    const std::vector< std::vector< std::vector<double> > > & 
//       image(wcsmap().image());
   const std::vector< std::vector<float> > & 
      solid_angles(wcsmap().solidAngles());
   const std::vector< std::vector< std::vector<float> > > & 
      image(wcsmap().image());
   // Compute exposures using exposure map energy grid
   std::vector<double> map_energies;
   expmap.getEnergies(map_energies);
   std::vector<double> map_exposures;

   for (int k(0); k < map_energies.size(); k++) {
      double my_exposure(0);
      for (size_t i(0); i < wcsmap().nxpix(); i++) {
         for (size_t j(0); j < wcsmap().nypix(); j++) {
            astro::SkyDir dir(wcsmap().skyDir(i+1, j+1));
            if (expmap.withinMapRadius(dir)) {
               my_exposure += (solid_angles[i][j]*image[0][j][i]
                               *expmap(dir, k));
            }
         }
      }
      map_exposures.push_back(my_exposure);
   }
   
   // Interpolate on the requested energy grid.
   for (size_t k(0); k < energies.size(); k++) {
      size_t indx;
      if (energies[k] <= map_energies.front()) {
         indx = 0;
      } else if (energies[k] >= map_energies.back()) {
         indx = map_energies.size() - 2;
      } else {
         indx = std::upper_bound(map_energies.begin(), map_energies.end(),
                                 energies[k]) - map_energies.begin() - 1;
      }
      exposure.push_back(interpolatePowerLaw(energies[k],
                                             map_energies[indx],
                                             map_energies[indx+1],
                                             map_exposures[indx],
                                             map_exposures[indx+1]));
   }
}

} // namespace Likelihood
