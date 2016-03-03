/**
 * @file MeanPsf.cxx
 * @brief Psf at a specific sky location averaged over an observation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/MeanPsf.cxx,v 1.34 2015/12/10 00:58:00 echarles Exp $
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
   computeExposure();
   if (s_separations.size() == 0) {
      createLogArray(1e-4, 70., 400, s_separations);
   }
   m_psfValues.reserve(m_energies.size()*s_separations.size());
   for (unsigned int k = 0; k < m_energies.size(); k++) {
      for (unsigned int j = 0; j < s_separations.size(); j++) {
         double value(0);
         std::map<unsigned int, irfInterface::Irfs *>::const_iterator 
            resp = m_observation.respFuncs().begin();
         for (; resp != m_observation.respFuncs().end(); ++resp) {
            int evtType = resp->second->irfID();
            Psf psf(s_separations[j], m_energies[k], evtType, m_observation);
            value += m_observation.expCube().value(m_srcDir, psf,
                                                   m_energies[k]);
         }
         if (m_exposure[k] > 0) {
            value /= m_exposure[k];
         } else {
            value = 0;
         }
         m_psfValues.push_back(value);
       }
   }
   // Ensure normalization and compute partial integrals.
   computePartialIntegrals();

   for (unsigned int k = 0; k < m_energies.size(); k++) {
     m_psfPeakValues.push_back((*this)(m_energies[k],0.0));
   }
}

void MeanPsf::computeExposure() {
   m_exposure.reserve(m_energies.size());
   for (size_t k(0); k < m_energies.size(); k++) {
      double value(0);
      std::map<unsigned int, irfInterface::Irfs *>::const_iterator
         resp = m_observation.respFuncs().begin();
      for (; resp != m_observation.respFuncs().end(); ++resp) {
         int evtType = resp->second->irfID();
         ExposureCube::Aeff aeff(m_energies[k], evtType, m_observation);
         value += m_observation.expCube().value(m_srcDir, aeff, m_energies[k]);
      }
      m_exposure.push_back(value);
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

double MeanPsf::peakValue(double energy) const {
   return st_facilities::Util::interpolate(m_energies, m_psfPeakValues, energy);
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

void MeanPsf::computePartialIntegrals() {
   m_partialIntegrals.clear();
   for (size_t k(0); k < m_energies.size(); k++) {
      std::vector<double> partialIntegral;
      partialIntegral.push_back(0);
      for (size_t j(0); j < s_separations.size()-1; j++) {
         size_t index(k*s_separations.size() + j);
         double theta1(s_separations.at(j)*M_PI/180.);
         double theta2(s_separations.at(j+1)*M_PI/180.);
         double y1(2.*M_PI*m_psfValues.at(index)*std::sin(theta1));
         double y2(2.*M_PI*m_psfValues.at(index+1)*std::sin(theta2));
         double slope((y2 - y1)/(theta2 - theta1));
         double intercept(y1 - theta1*slope);
         double value(slope*(theta2*theta2 - theta1*theta1)/2.
                      + intercept*(theta2 - theta1));
         partialIntegral.push_back(partialIntegral.back() + value);
      }
      // Ensure normalizations of differential and integral arrays.
      for (size_t j(0); j < s_separations.size(); j++) {
         size_t index(k*s_separations.size() + j);
         m_psfValues.at(index) /= partialIntegral.back();
         partialIntegral[j] /= partialIntegral.back();
      }
      m_partialIntegrals.push_back(partialIntegral);
   }
}

double MeanPsf::integral(double angle, double energy) const {
   size_t k(std::upper_bound(m_energies.begin(), m_energies.end(),
                             energy) - m_energies.begin() - 1);
   if (k < 0 || k > static_cast<int>(m_energies.size()-1)) {
      std::ostringstream what;
      what << "MeanPsf::integral: energy " << energy << " out-of-range.";
      throw std::out_of_range(what.str());
   }

   if (angle < s_separations.front()) {
      return 0;
   } else if (angle >= s_separations.back()) {
      return 1;
   }
   size_t j(std::upper_bound(s_separations.begin(), s_separations.end(),
                             angle) - s_separations.begin() - 1);

   size_t index(k*s_separations.size() + j);
   double theta1(s_separations[j]*M_PI/180.);
   double theta2(s_separations[j+1]*M_PI/180.);
   double y1(2.*M_PI*m_psfValues.at(index)*std::sin(theta1));
   double y2(2.*M_PI*m_psfValues.at(index+1)*std::sin(theta2));
   double slope((y2 - y1)/(theta2 - theta1));
   double intercept(y1 - theta1*slope);
   double theta(angle*M_PI/180.);
   double value(slope*(theta*theta - theta1*theta1)/2. 
                + intercept*(theta - theta1));
   double integral1(m_partialIntegrals.at(k).at(j) + value);

   if (k == m_energies.size()-1) {
      // energy is at upper boundary of m_energies vector so just
      // use the psf integral evaluated at that bound.
      return integral1;
   }
   index += s_separations.size();
   y1 = 2.*M_PI*m_psfValues.at(index)*std::sin(theta1);
   y2 = 2.*M_PI*m_psfValues.at(index+1)*std::sin(theta2);
   slope = (y2 - y1)/(theta2 - theta1);
   intercept = y1 - theta1*slope;
   value = slope*(theta*theta - theta1*theta1)/2. + intercept*(theta - theta1);
   double integral2(m_partialIntegrals.at(k+1).at(j) + value);
   double final_value((energy - m_energies.at(k))
                      /(m_energies.at(k+1) - m_energies.at(k))
                      *(integral2 - integral1) + integral1);
   return final_value;
}

// double MeanPsf::containmentRadius(double energy, double frac) const {
//    if (energy < m_energies.front() || energy > m_energies.back()) {
//       std::ostringstream message;
//       message << "MeanPsf::containmentRadius: selected energy, "
//               << energy << " is out-of-range.";
//       throw std::out_of_range(message.str());
//    }
//    size_t k(std::upper_bound(m_energies.begin(), m_energies.end(), energy)
//             - m_energies.begin() - 1);
//    if (frac <= 0) {
//       return 0;
//    } else if (frac >= 1) {
//       return 1;
//    }
//    size_t j1(std::upper_bound(m_partialIntegrals[k].begin(),
//                               m_partialIntegrals[k].end(), frac)
//              - m_partialIntegrals[k].begin() - 1);
//    size_t j2(std::upper_bound(m_partialIntegrals[k+1].begin(),
//                               m_partialIntegrals[k+1].end(), frac)
//              - m_partialIntegrals[k+1].begin() - 1);
//    double angle1( (frac - m_partialIntegrals[k][j1])
//                   /(m_partialIntegrals[k][j1+1] - m_partialIntegrals[k][j1])
//                   *(s_separations[j1+1] - s_separations[j1])
//                   + s_separations[j1] );
//    double angle2( (frac - m_partialIntegrals[k+1][j2])
//                   /(m_partialIntegrals[k+1][j2+1] - m_partialIntegrals[k+1][j2])
//                   *(s_separations[j2+1] - s_separations[j2])
//                   + s_separations[j2] );
//    double value( (energy - m_energies[k])/(m_energies[k+1] - m_energies[k])
//                  *(angle2 - angle1) + angle1);
//    return value;
// }

double MeanPsf::derivative(double angle, double energy) const {
   int j = std::upper_bound(s_separations.begin(), s_separations.end(),
                            angle) - s_separations.begin() - 1;
   int k = std::upper_bound(m_energies.begin(), m_energies.end(),
                            energy) - m_energies.begin() - 1;
   std::ostringstream what;
   if (j < 0 || j >= static_cast<int>(s_separations.size()-1)) {
      what << "MeanPsf::derivative: angle " << angle << " out-of-range.";
      throw std::out_of_range(what.str());
   }
   if (k < 0 || k > static_cast<int>(m_energies.size()-1)) {
      what << "MeanPsf::derivative: energy " << energy << " out-of-range.";
      throw std::out_of_range(what.str());
   }
   int index = k*s_separations.size() + j;
   double value = ((m_psfValues.at(index+1) - m_psfValues.at(index))
                   /(s_separations[j+1] - s_separations[j]));
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

void MeanPsf::getImage(double energy, double lon0, double lat0,
		       const std::vector<std::pair<double,double> >& dirs,
		       std::vector<double> & image) const {
    astro::SkyDir center(lon0, lat0);
    image.clear();
    image.reserve(dirs.size());
    for (unsigned int i = 0; i < dirs.size(); i++) {
        const std::pair<double,double>& thePair = dirs[i];
        astro::SkyDir my_dir(thePair.first,thePair.second);
	double dist = my_dir.difference(center)*180./M_PI;
	image.push_back(this->operator()(energy, dist, 0));
    }
}



void MeanPsf::createLogArray(double xmin, double xmax, unsigned int npts,
                             std::vector<double> &xx) const {
   xx.clear();
   xx.reserve(npts+1);
   xx.push_back(0);
   double xstep = log(xmax/xmin)/(npts - 1.);
   for (unsigned int i = 0; i < npts; i++) {
      xx.push_back(xmin*exp(i*xstep));
   }
}

double MeanPsf::Psf::value(double cosTheta, double phi) const {
   double inclination = acos(cosTheta)*180./M_PI;
   double epoch((m_observation.expCube().tstart() 
                 + m_observation.expCube().tstop())/2.);
   std::map<unsigned int, irfInterface::Irfs *>::const_iterator respIt 
      = m_observation.respFuncs().begin();
   for ( ; respIt != m_observation.respFuncs().end(); ++respIt) {
      if (respIt->second->irfID() == m_evtType) {  
         irfInterface::IAeff * aeff = respIt->second->aeff();
         irfInterface::IPsf * psf = respIt->second->psf();
         double aeffValue = aeff->value(m_energy, inclination, phi, epoch);
         double psfValue = psf->value(m_separation, m_energy, inclination, 
                                      phi, epoch);
         double psf_val = aeffValue*psfValue;
         if (psf_val < 0) {
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
