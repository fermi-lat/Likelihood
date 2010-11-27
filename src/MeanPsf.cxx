/**
 * @file MeanPsf.cxx
 * @brief Psf at a specific sky location averaged over an observation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/MeanPsf.cxx,v 1.23 2010/11/25 12:52:38 cohen Exp $
 */

#include <cmath>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "st_stream/StreamFormatter.h"

#include "st_facilities/Util.h"

#include "Likelihood/MeanPsf.h"

namespace Likelihood {

std::vector<double> MeanPsf::s_separations;

void MeanPsf::init() {
   if (s_separations.size() == 0) {
      createLogArray(1e-4, 70., 200, s_separations);
   }

   m_psfValues.reserve(m_energies.size()*s_separations.size());
   m_exposure.reserve(m_energies.size());
   for (unsigned int k = 0; k < m_energies.size(); k++) {
      for (unsigned int j = 0; j < s_separations.size(); j++) {
         double expsr_val(0);
         double psf_val(0);
         std::map<unsigned int, irfInterface::Irfs *>::const_iterator 
            resp = m_observation.respFuncs().begin();
         for (; resp != m_observation.respFuncs().end(); ++resp) {
            int evtType = resp->second->irfID();
//            Aeff aeff(m_energies[k], evtType, m_observation);
            ExposureCube::Aeff aeff(m_energies[k], evtType, m_observation);
            expsr_val += m_observation.expCube().value(m_srcDir, aeff,
                                                       m_energies[k]);
            Psf psf(s_separations[j], m_energies[k], evtType, m_observation);
            psf_val += m_observation.expCube().value(m_srcDir, psf,
                                                     m_energies[k]);
         }
         if (j == 0) {
            m_exposure.push_back(expsr_val);
         }
         if (expsr_val > 0) {
            psf_val /= expsr_val;
         } else {
            psf_val = 0;
         }
         m_psfValues.push_back(psf_val);
       }
   }
}

double MeanPsf::operator()(double energy, double theta, double phi) const {
   (void)(phi);
   if (m_energies.size() == 0) {
      throw std::runtime_error("MeanPsf::operator(): Cannot call a " + 
                               std::string("default object."));
   }
   if (energy < m_energies.front() || energy > m_energies.back()) {
      return 0;
   }
   if (theta > s_separations.back()) {
      return 0;
   }
   return st_facilities::Util::bilinear(m_energies, energy, 
                                        s_separations, theta, 
                                        m_psfValues);
}

void MeanPsf::write(const std::string & filename) const {
   std::ofstream output_file(filename.c_str());
   std::vector<double>::const_iterator sep = s_separations.begin();
   for (unsigned int j = 0; sep != s_separations.end(); ++sep, j++) {
      output_file << *sep << "  ";
      for (unsigned int k = 0; k < m_energies.size(); k++) {
         unsigned int indx = k*s_separations.size() + j;
         output_file << m_psfValues.at(indx) << "  ";
      }
      output_file << std::endl;
   }
   output_file.close();
}

double MeanPsf::exposure(double energy) const {
   return st_facilities::Util::interpolate(m_energies, m_exposure, energy);
}

double MeanPsf::integral(double angle, double energy) const {
   int sepmax = std::upper_bound(s_separations.begin(), s_separations.end(),
                                 angle) - s_separations.begin();
   int k = std::upper_bound(m_energies.begin(), m_energies.end(),
                            energy) - m_energies.begin() - 1;
   if (sepmax <= 0) {
      return 0;
   } else if (sepmax > static_cast<int>(s_separations.size()-1)) {
      return 1;
   }
   if (k < 0 || k > static_cast<int>(m_energies.size()-1)) {
      std::ostringstream what;
      what << "MeanPsf::integral: energy " << energy << " out-of-range.";
      throw std::out_of_range(what.str());
   }
   double theta1;
   double theta2;
   double integral1(0);
   double integral2(0);
   for (int j = 0; j < sepmax; j++) {
      int indx1 = k*s_separations.size() + j;
      theta1 = s_separations.at(j)*M_PI/180.;
      theta2 = s_separations.at(j+1)*M_PI/180.;
      integral1 += ((m_psfValues.at(indx1)*std::sin(theta1)
                     + m_psfValues.at(indx1+1)*std::sin(theta2))/2.
                    *(theta2 - theta1));

      int indx2 = indx1 + s_separations.size();
      integral2 += ((m_psfValues.at(indx2)*std::sin(theta1)
                     + m_psfValues.at(indx2+1)*std::sin(theta2))/2.
                    *(theta2 - theta1));
   }
   theta1 = s_separations.at(sepmax-1)*M_PI/180.;
   theta2 = angle*M_PI/180.;
   integral1 += (operator()(m_energies.at(k), angle)*std::sin(theta2) +
                 operator()(m_energies.at(k), s_separations.at(sepmax-1))
                 *std::sin(theta1))/2.*(theta2 - theta1);

   integral2 += (operator()(m_energies.at(k+1), angle)*std::sin(theta2) +
                 operator()(m_energies.at(k+1), s_separations.at(sepmax-1))
                 *std::sin(theta1))/2.*(theta2 - theta1);

   double value = 2.*M_PI*((energy - m_energies.at(k))/
                           (m_energies.at(k+1) - m_energies.at(k))
                           *(integral2 - integral1) + integral1);
   return value;
}

void MeanPsf::getImage(double energy, double lon0, double lat0,
                       const std::vector<double> lons,
                       const std::vector<double> lats,
                       std::vector< std::vector<double> > & image) const {
   astro::SkyDir center(lon0, lat0);
   image.clear();
   image.reserve(lons.size());
   for (unsigned int i = 0; i < lons.size(); i++) {
      std::vector<double> row;
      for (unsigned int j = 0; j < lats.size(); j++) {
         astro::SkyDir my_dir(lons.at(i), lats.at(j));
         double dist = my_dir.difference(center)*180./M_PI;
         row.push_back(this->operator()(energy, dist, 0));
      }
      image.push_back(row);
   }
}

void MeanPsf::createLogArray(double xmin, double xmax, unsigned int npts,
                             std::vector<double> &xx) const {
   xx.clear();
   xx.reserve(npts);
   double xstep = log(xmax/xmin)/(npts - 1.);
   for (unsigned int i = 0; i < npts; i++) {
      xx.push_back(xmin*exp(i*xstep));
   }
}

double MeanPsf::Psf::operator()(double cosTheta, double phi) const {
   double inclination = acos(cosTheta)*180./M_PI;
   if (inclination > 70.) {
      return 0;
   }
   std::map<unsigned int, irfInterface::Irfs *>::const_iterator respIt 
      = m_observation.respFuncs().begin();
   for ( ; respIt != m_observation.respFuncs().end(); ++respIt) {
      if (respIt->second->irfID() == m_evtType) {  
         irfInterface::IAeff * aeff = respIt->second->aeff();
         irfInterface::IPsf * psf = respIt->second->psf();
         double aeffValue = aeff->value(m_energy, inclination, phi);
         double psfValue = psf->value(m_separation, m_energy, inclination, phi);
         double psf_val = aeffValue*psfValue;
         if (psf_val < 0) {
            if (inclination > 69.) {  // ugly kluge
               return 0;
            }
            st_stream::StreamFormatter formatter("MeanPsf", "operator()", 4);
            formatter.info() << "separation: " << m_separation << "  "
                             << "energy: " << m_energy << "  "
                             << "inclination: " <<inclination << "  "
                             << "phi: " << phi << std::endl;
            throw std::runtime_error("MeanPsf::Psf::operator(): psf_val < 0");
         }
         return psf_val;
      }
   }
   return 0;
}

} // namespace Likelihood
