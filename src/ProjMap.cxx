/**
 * @file WcsMap2.cxx
 * @brief A map with reference point centered on the image and that
 * uses WCS projections for indexing its internal representation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/echarles/healpix_changes/Likelihood/src/ProjMap.cxx,v 1.3 2015/12/02 00:53:06 echarles Exp $
 */

#include <cmath>

#include <algorithm>
#include <iostream>
#include <stdexcept>

#include "st_stream/StreamFormatter.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Image.h"

#include "st_facilities/Util.h"

#include "astro/SkyProj.h"
#include "astro/ProjBase.h"

#include "Likelihood/SkyDirArg.h"
#include "Likelihood/ProjMap.h"

namespace Likelihood {

ProjMap::ProjMap() 
  : m_filename(),
    m_refDir(0., 0.),
    m_mapRadius(0.),
    m_proj(0),
    m_interpolate(true),
    m_mapIntegral(0),
    m_enforceEnergyRange(false),
    m_extrapolated(0){
}


ProjMap::ProjMap(const std::string & filename,
                 bool interpolate,
                 bool enforceEnergyRange) 
  : m_filename(filename),
    m_refDir(0.,0.),
    m_mapRadius(0.),
    m_proj(0),
    m_interpolate(interpolate),
    m_mapIntegral(0),
    m_enforceEnergyRange(enforceEnergyRange),
    m_extrapolated(0) {
}
  


ProjMap::~ProjMap() {
   delete m_proj;
   st_stream::StreamFormatter formatter("ProjMap", "", 2);
   if (m_extrapolated > 0) {
      formatter.info(3) << "ProjMap: extrapolated beyond the maximum "
                        << "energy for map cube file " << m_filename << " "
                        << m_extrapolated << " times.\n";
   }
}

ProjMap::ProjMap(const ProjMap & rhs) 
  :  m_filename(rhs.m_filename),
     m_refDir(rhs.m_refDir),
     m_mapRadius(rhs.m_mapRadius),
     m_proj(rhs.m_proj ? rhs.m_proj->clone() : 0),
     m_energies(rhs.m_energies),
     m_interpolate(rhs.m_interpolate),
     m_enforceEnergyRange(rhs.m_enforceEnergyRange),
     m_extrapolated(rhs.m_extrapolated), 
     m_mapIntegral(rhs.m_mapIntegral),
     m_mapIntegrals(rhs.m_mapIntegrals){
}

ProjMap::ProjMap(const ProjMap & rhs, const double & energy) 
  :  m_filename(rhs.m_filename),
     m_refDir(rhs.m_refDir),
     m_mapRadius(rhs.m_mapRadius),
     m_proj(rhs.m_proj ? rhs.m_proj->clone() : 0),
     m_energies(1,energy),
     m_interpolate(rhs.m_interpolate),
     m_enforceEnergyRange(rhs.m_enforceEnergyRange),
     m_extrapolated(rhs.m_extrapolated), 
     m_mapIntegral(0.),
     m_mapIntegrals(1,0.){
}

ProjMap & ProjMap::operator=(const ProjMap & rhs) {
   if (this != &rhs) {
      m_filename = rhs.m_filename;
      m_refDir = rhs.m_refDir;
      m_mapRadius = rhs.m_mapRadius;
      delete m_proj;
      m_proj = rhs.m_proj ? rhs.m_proj->clone() : 0;
      m_energies = rhs.m_energies;
      m_interpolate = rhs.m_interpolate;
      m_enforceEnergyRange = rhs.m_enforceEnergyRange;
      m_mapIntegral = rhs.m_mapIntegral;
      m_mapIntegrals = rhs.m_mapIntegrals;
      m_extrapolated = rhs.m_extrapolated;
   }
   return *this;
}


astro::SkyDir ProjMap::skyDir(double ilon, double ilat) const {
   if ( m_proj == 0 ) {
     throw std::runtime_error("ProjMap: Projection is not initialized");
   }
   // EAC -> debugging
   std::pair<double, double> coords(0.,0.);
   static int count(0);
   try {
     coords = m_proj->pix2sph(ilon, ilat);
   } catch (...) {
     if ( count % 1000 == 0 ) {
       std::cerr << "ProjMap::skyDir " << ilon << ' ' << ilat << std::endl;       
     }
     count++;
   }

   return astro::SkyDir(coords.first, coords.second, 
			m_proj->isGalactic() ? astro::SkyDir::GALACTIC : astro::SkyDir::EQUATORIAL );
}

bool ProjMap::withinMapRadius(const astro::SkyDir & dir) const {
   if (dir.difference(m_refDir) <= m_mapRadius) {
      return true;
   }
   return false;
}



double ProjMap::mapIntegral(double energy) const {
   check_energy(energy);

   int k(0);
   if (m_energies.size() > 1) {
      k = std::upper_bound(m_energies.begin(), m_energies.end(), energy)
         - m_energies.begin() - 1;
   }
   if (energy == m_energies.at(k)) {
      return m_mapIntegrals.at(k);
   }
   if (k == m_energies.size() - 1) {
      k = m_energies.size() - 2;
   }

   double value = (m_mapIntegrals[k]*
                   std::exp((std::log(energy/m_energies[k]))
                            /(std::log(m_energies[k+1]/m_energies[k]))
                            *std::log(m_mapIntegrals[k+1]/m_mapIntegrals[k])));
   return value;
}

void ProjMap::setProjInfo(const astro::SkyDir& dir, const astro::ProjBase& proj) {
  // Take ownershipe of the projection
  m_proj = const_cast<astro::ProjBase*>(&proj);
  m_refDir = dir;
}

void ProjMap::check_energy_index(int k) const {
   if ( k < 0 || k >  m_energies.size() ) {
      throw std::runtime_error("ProjMap: Requested energy index is "
                               "out-of-range.");
   }
}

void ProjMap::check_energy(double energy) const {
   if (m_enforceEnergyRange && 
       (energy < m_energies.front() || energy > m_energies.back())) {
      throw std::runtime_error("ProjMap: Requested energy is out-of-range.");
   }
}

double ProjMap::interpolatePowerLaw(double x, double x1, double x2,
                                    double y1, double y2) {
   if (y1 == 0 && y2 == 0) {
      return 0;
   }
   if (y1 == 0 || y2 == 0) {
      if ( (x < x1 && y1 == 0) || (x > x2 && y2 == 0) ) {
         // Linear extrapolation in these cases would produce
         // negative values.
         std::ostringstream message;
         message << "ProjMap::interpolatePowerLaw: "
                 << "linear extrapolation selected for "
                 << "zero-valued ordinates.\n"
                 << "x = " << x << ", "
                 << "x1 = " << x1 << ", "
                 << "x2 = " << x2 << ", "
                 << "y1 = " << y1 << ", "
                 << "y2 = " << y2 << std::endl;
         throw std::runtime_error(message.str());
      }
      // Use linear interpolation
      double value((x - x1)/(x2 - x1)*(y2 - y1) + y1);
      return value;
   }
   if (x1 <= 0 || x2 <= 0 || y1 <= 0 || y2 <= 0) {
      std::ostringstream message;
      message << "ProjMap::interpolatePowerLaw:\n"
              << "abscissa or ordinate values found that are <= 0: "
              << "x1 = " << x1 << ", "
              << "x2 = " << x2 << ", "
              << "y1 = " << y1 << ", "
              << "y2 = " << y2 << std::endl;
      throw std::runtime_error(message.str());
   }
   double gamma = std::log(y2/y1)/std::log(x2/x1);
   // double n0 = y1/std::pow(x1, gamma);
   // return n0*std::pow(x, gamma);
   return y1*std::pow(x/x1, gamma);
}

} // namespace Likelihood
