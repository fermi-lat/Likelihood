/**
 * @file MeanPsf.cxx
 * @brief Psf averaged over an observation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/MeanPsf.cxx,v 1.9 2004/11/28 06:58:21 jchiang Exp $
 */

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "st_facilities/Util.h"

#include "Likelihood/ExposureCube.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/ResponseFunctions.h"

#include "Verbosity.h"

namespace Likelihood {

std::vector<double> MeanPsf::s_separations;

void MeanPsf::init() {
   if (ExposureCube::instance() == 0) {
      throw std::runtime_error("Likelihood::MeanPsf: No exposure hypercube "
                               + std::string("file."));
   }
   if (s_separations.size() == 0) {
      createLogArray(1e-4, 70., 200, s_separations);
   }

   m_psfValues.reserve(m_energies.size()*s_separations.size());
   for (unsigned int k = 0; k < m_energies.size(); k++) {
      for (unsigned int j = 0; j < s_separations.size(); j++) {
         double aeff_val(0);
         double psf_val(0);
         for (int evtType = 0; evtType < 2; evtType++) {
            Aeff aeff(m_energies[k], evtType);
            aeff_val += ExposureCube::instance()->value(m_srcDir, aeff);;
            Psf psf(s_separations[j], m_energies[k], evtType);
            psf_val += ExposureCube::instance()->value(m_srcDir, psf);
         }
         if (aeff_val > 0) {
            psf_val /= aeff_val;
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

void MeanPsf::createLogArray(double xmin, double xmax, unsigned int npts,
                             std::vector<double> &xx) const {
   xx.clear();
   xx.reserve(npts);
   double xstep = log(xmax/xmin)/(npts - 1.);
   for (unsigned int i = 0; i < npts; i++) {
      xx.push_back(xmin*exp(i*xstep));
   }
}

double MeanPsf::Psf::s_phi(0);

double MeanPsf::Psf::operator()(double cosTheta) const {
   double inclination = acos(cosTheta)*180./M_PI;
   if (inclination > 70.) {
      return 0;
   }
   std::map<unsigned int, irfInterface::Irfs *>::const_iterator respIt 
      = ResponseFunctions::instance()->begin();
   for ( ; respIt != ResponseFunctions::instance()->end(); ++respIt) {
      if (respIt->second->irfID() == m_evtType) {  
         irfInterface::IAeff * aeff = respIt->second->aeff();
         irfInterface::IPsf * psf = respIt->second->psf();
         double aeffValue = aeff->value(m_energy, inclination, s_phi);
         double psfValue = psf->value(m_separation, m_energy, inclination, 
                                      s_phi);
         double psf_val = aeffValue*psfValue;
         if (psf_val < 0) {
            if (inclination > 69.) {  // ugly kluge
               return 0;
            }
            if (print_output(4)) {
               std::cerr << "separation: " << m_separation << "  "
                         << "energy: " << m_energy << "  "
                         << "inclination: " <<inclination << "  "
                         << "phi: " << s_phi << std::endl;
            }
            throw std::runtime_error("MeanPsf::Psf::operator(): psf_val < 0");
         }
         return psf_val;
      }
   }
   return 0;
}

double MeanPsf::Aeff::s_phi(0);

double MeanPsf::Aeff::operator()(double cosTheta) const {
   double inclination = acos(cosTheta)*180./M_PI;
   std::map<unsigned int, irfInterface::Irfs *>::const_iterator respIt 
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
