/**
 * @file MapCubeFunction2.cxx
 * @brief Encapsulation of 3D FITS image of a diffuse source with 
 * position-dependent spectral variation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/echarles/healpix_changes/Likelihood/src/MapCubeFunction2.cxx,v 1.5 2015/11/25 18:52:42 echarles Exp $
 */

#include <cmath>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "facilities/Util.h"

#include "st_stream/StreamFormatter.h"

#include "Likelihood/ExposureMap.h"
#include "Likelihood/MapCubeFunction2.h"
#include "Likelihood/SkyDirArg.h"
#include "Likelihood/WcsMap2.h"
#include "Likelihood/HealpixProjMap.h"

namespace Likelihood {

MapCubeFunction2::MapCubeFunction2() 
   : optimizers::Function("MapCubeFunction", 1, "Normalization"), MapBase() {
  addParam("Normalization", 1, false);

}

MapCubeFunction2::MapCubeFunction2(const std::string & fitsFile) 
   : optimizers::Function("MapCubeFunction", 1, "Normalization"), 
     MapBase(fitsFile) {
  addParam("Normalization", 1, false);
}

MapCubeFunction2::MapCubeFunction2(const MapCubeFunction2 & rhs)
   : optimizers::Function(rhs), MapBase(rhs) {}

MapCubeFunction2 & MapCubeFunction2::operator=(const MapCubeFunction2 & rhs) {
   if (this != &rhs) {
      optimizers::Function::operator=(rhs);
      MapBase::operator=(rhs);
   }
   return *this;
}

MapCubeFunction2::~MapCubeFunction2() {
}

double MapCubeFunction2::value(const optimizers::Arg & xarg) const {
   const SkyDirArg & dir = dynamic_cast<const SkyDirArg &>(xarg);
   double energy = dir.energy();
   double value = projmap().operator()(dir(), energy);
   return value*getParam("Normalization").getTrueValue();
}


double MapCubeFunction2::mapIntegral() const {
   return projmap().mapIntegral();
}

double MapCubeFunction2::mapIntegral(double energy) const {
   return projmap().mapIntegral(energy);
}

void MapCubeFunction2::
integrateSpatialDist(const std::vector<double> & energies,
		     const ExposureMap & expmap,
		     std::vector<double> & exposure) const {
   // EAC, switch based on projection type
   const ProjMap& projMap = projmap();
   switch ( projMap.getProj()->method() ) {
   case astro::ProjBase::WCS:
     return integrateSpatialDist_wcs(energies,expmap,static_cast<const WcsMap2&>(projMap),exposure);
   case astro::ProjBase::HEALPIX:
     return integrateSpatialDist_healpix(energies,expmap,static_cast<const HealpixProjMap&>(projMap),exposure);
   default:
     break;
   }
   std::string errMsg("Unrecognized projection type for MapCubeFunction2: ");
   errMsg += fitsFile();
   throw std::runtime_error(errMsg);
   return;
}


void MapCubeFunction2::
integrateSpatialDist_healpix(const std::vector<double> & energies,
			     const ExposureMap & expmap,
			     const HealpixProjMap& wcsmap,
			     std::vector<double> & exposure) const {
   exposure.clear();

   const double& solidAngle = wcsmap.solidAngleHealpix();
   const std::vector< Healpix_Map<float>  >& image =  wcsmap.image();

   // Compute exposures using exposure map energy grid
   std::vector<double> map_energies;
   expmap.getEnergies(map_energies);
   std::vector<double> map_exposures;

   for (int k(0); k < map_energies.size(); k++) {
      double my_exposure(0);
      for (size_t i(0); i < wcsmap.nPixels(); i++) {
         astro::SkyDir dir(wcsmap.skyDir(i, 0));
	 if (expmap.withinMapRadius(dir)) {
	     my_exposure += (solidAngle
			     *wcsmap(dir, map_energies[k])
			     *expmap(dir, k));
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

  return;
}

void MapCubeFunction2::
integrateSpatialDist_wcs(const std::vector<double> & energies,
			 const ExposureMap & expmap,
			 const WcsMap2& wcsmap,
			 std::vector<double> & exposure) const {
   exposure.clear();
//    const std::vector< std::vector<double> > & 
//       solid_angles(projmap().solidAngles());
//    const std::vector< std::vector< std::vector<double> > > & 
//       image(projmap().image());
   const std::vector< std::vector<float> > & 
      solid_angles(wcsmap.solidAngles());
   const std::vector< std::vector< std::vector<float> > > & 
      image(wcsmap.image());

   // Compute exposures using exposure map energy grid
   std::vector<double> map_energies;
   expmap.getEnergies(map_energies);
   std::vector<double> map_exposures;

   for (int k(0); k < map_energies.size(); k++) {
      double my_exposure(0);
      for (size_t i(0); i < wcsmap.nxpix(); i++) {
         for (size_t j(0); j < wcsmap.nypix(); j++) {
            astro::SkyDir dir(wcsmap.skyDir(i+1, j+1));
            if (expmap.withinMapRadius(dir)) {
               my_exposure += (solid_angles[i][j]
                               *wcsmap(dir, map_energies[k])
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
