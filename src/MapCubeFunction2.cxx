/**
 * @file MapCubeFunction2.cxx
 * @brief Encapsulation of 3D FITS image of a diffuse source with 
 * position-dependent spectral variation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/MapCubeFunction2.cxx,v 1.5 2012/01/06 07:11:59 jchiang Exp $
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

namespace Likelihood {

MapCubeFunction2::MapCubeFunction2() 
   : optimizers::Function(), MapBase() {
   init();
}

MapCubeFunction2::MapCubeFunction2(const std::string & fitsFile) 
   : optimizers::Function(), MapBase(fitsFile) {
   init();
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

double MapCubeFunction2::value(optimizers::Arg & xarg) const {
   SkyDirArg & dir = dynamic_cast<SkyDirArg &>(xarg);
   double energy = dir.energy();
   double value = wcsmap().operator()(dir(), energy);
   return value*getParam("Normalization").getTrueValue();
}

void MapCubeFunction2::init() {
   setMaxNumParams(1);
// Leave this parameter fixed, modifying the overall normalization
// via a ConstantValue function as the spectral component.
   addParam("Normalization", 1, false);

   m_funcType = Addend;
   m_argType = "";
   m_genericName = "MapCubeFunction";
   m_normParName = "Normalization";
}

double MapCubeFunction2::mapIntegral() const {
   return wcsmap().mapIntegral();
}

double MapCubeFunction2::mapIntegral(double energy) const {
   return wcsmap().mapIntegral(energy);
}

void MapCubeFunction2::
integrateSpatialDist(const std::vector<double> & energies,
                     const ExposureMap & expmap,
                     std::vector<double> & exposure) const {
   exposure.clear();
   const std::vector< std::vector<double> > & 
      solid_angles(wcsmap().solidAngles());
   const std::vector< std::vector< std::vector<double> > > & 
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
               my_exposure += (solid_angles[i][j]
                               *wcsmap()(dir, map_energies[k])
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
