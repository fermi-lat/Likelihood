/**
 * @file BinnedExposure.cxx
 * @brief Integral of effective area over time for the entire sky at
 * various energies.
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>
#include <algorithm>
#include <stdexcept>

#include "Likelihood/BinnedExposure.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/ExposureCube.h"
#include "Likelihood/ResponseFunctions.h"

namespace Likelihood {

BinnedExposure::BinnedExposure(const std::vector<double> & energies) 
   : m_energies(energies) {
   computeMap();
}

BinnedExposure::BinnedExposure(const std::string & filename) {
   CountsMap exposureMap(filename);
   exposureMap.getAxisVector(0, m_ras);
   exposureMap.getAxisVector(1, m_decs);
   exposureMap.getAxisVector(2, m_energies);
   m_exposureMap = exposureMap.data();
}

double BinnedExposure::operator()(double energy, double ra, double dec) const {
   std::vector<double>::const_iterator ie = std::find(m_energies.begin(),
                                                      m_energies.end(),
                                                      energy);
   if (ie == m_energies.end()) {
      throw std::runtime_error("BinnedExposure::operator(): The energy " +
                               std::string("selected is not available."));
   }
   unsigned int k = ie - m_energies.begin();
   if (ra < 0) ra += 360.; 
   if (ra > 360.) ra = fmod(ra, 360.);
   unsigned int i = findIndex(m_ras.begin(), m_ras.end(), ra);
   unsigned int j = findIndex(m_decs.begin(), m_decs.end(), dec);
   unsigned int indx = k*m_ras.size()*m_decs.size() + j*m_ras.size() + i;
   return m_exposureMap.at(indx);
}

unsigned int 
BinnedExposure::findIndex(std::vector<double>::const_iterator begin,
                          std::vector<double>::const_iterator end,
                          double value) const {
   std::vector<double>::const_iterator it 
      = std::upper_bound(begin, end, value);
   if (it == end) {
      throw std::range_error("BinnedExposure::findIndex: out of range value.");
   }
   int indx = it - begin - 1;
   return indx;
}

void BinnedExposure::computeMap() {
   linearArray(0., 360., 361, m_ras);
   linearArray(-90., 90., 181, m_decs);

   m_exposureMap.resize(m_ras.size()*m_decs.size()*m_energies.size(), 0);
   int iter(0);
   for (unsigned int j = 0; j < m_decs.size(); j++) {
      for (unsigned int i = 0; i < m_ras.size(); i++) {
         if ( (iter % ((m_decs.size()*m_ras.size())/20)) == 0 ) {
            std::cerr << ".";
         }
         astro::SkyDir dir(m_ras[i], m_decs[j], astro::SkyDir::EQUATORIAL);
         for (unsigned int k = 0; k < m_energies.size(); k++) {
            unsigned int indx 
               = k*m_ras.size()*m_decs.size() + j*m_ras.size() + i;
            for (int evtType = 0; evtType < 2; evtType++) {
               Aeff aeff(m_energies[k], evtType);
               m_exposureMap.at(indx)
                  += ExposureCube::instance()->value(dir, aeff);
            }
         }
         iter++;
      }
   }
}

void BinnedExposure::linearArray(double xmin, double xmax, unsigned int npts,
                                 std::vector<double> &xx) const {
   double xstep = (xmax - xmin)/(npts - 1.);
   for (unsigned int i = 0; i < npts; i++) {
      xx.push_back(i*xstep + xmin);
   }
}

double BinnedExposure::Aeff::s_phi(0);

double BinnedExposure::Aeff::operator()(double cosTheta) const {
   double inclination = acos(cosTheta);
   std::map<unsigned int, irfInterface::Irfs *>::iterator respIt 
      = ResponseFunctions::instance()->begin();
   for ( ; respIt != ResponseFunctions::instance()->end(); ++respIt) {
      if (respIt->second->irfID() == m_evtType) {  
         irfInterface::IAeff * aeff = respIt->second->aeff();
         double aeff_val = aeff->value(m_energy, inclination, s_phi);
         return aeff_val;
      }
   }
   return 0;
}

} // namespace Likelihood
