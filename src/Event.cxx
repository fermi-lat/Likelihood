/** 
 * @file Event.cxx
 * @brief Event class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Event.cxx,v 1.70 2009/02/23 00:38:11 jchiang Exp $
 */

#include <cctype>

#include <algorithm>
#include <deque>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <utility>

#include "Likelihood/DiffuseSource.h"
#include "Likelihood/Event.h"
#include "Likelihood/EquinoxRotation.h"
#include "Likelihood/Exception.h"
#include "Likelihood/MapBase.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/ScData.h"
#include "Likelihood/SkyDirArg.h"
#include "Likelihood/Source.h"
#include "Likelihood/TrapQuad.h"

#include "st_facilities/GaussianQuadrature.h"
#include "Likelihood/DiffRespIntegrand.h"

#include "LogNormalMuDist.h"

namespace {
   double my_acos(double mu) {
      if (mu > 1) {
         return 0;
      } else if (mu < -1) {
         return M_PI;
      } else {
         return acos(mu);
      }
   }
}

namespace Likelihood {

std::vector<double> Event::s_mu;
std::vector<double> Event::s_mu_2;
std::vector<double> Event::s_phi;
bool Event::s_haveSourceRegionData(false);

Event::Event() : m_respName(""), m_modelSum(0), m_ctbclasslevel(0) {}

Event::Event(double ra, double dec, double energy, double time, 
             const astro::SkyDir & scZAxis, const astro::SkyDir & scXAxis, 
             double muZenith, bool useEdisp, const std::string & respName,
             int type) 
   : m_appDir(astro::SkyDir(ra, dec)), m_energy(energy), m_arrTime(time),
     m_muZenith(muZenith), m_type(type), m_scDir(scZAxis), m_scXDir(scXAxis),
     m_useEdisp(useEdisp), m_respName(respName), m_modelSum(0),
     m_ctbclasslevel(0) {
   if (m_useEdisp) {
// For <15% energy resolution, consider true energies over the range
// (0.55, 1.45)*m_energy, i.e., nominally a >3-sigma range about the
// apparent energy.
//       int npts(100);
      int npts(30);
      double emin = 0.55*m_energy;
      double emax = 1.45*m_energy;
      m_estep = (emax - emin)/(npts-1.);
      m_trueEnergies.reserve(npts);
      for (int i = 0; i < npts; i++) {
         m_trueEnergies.push_back(m_estep*i + emin);
      }
   } else {
// To mimic infinite energy resolution, we create a single element
// vector containing the apparent energy.
      m_trueEnergies.push_back(m_energy);
   }
}

double Event::diffuseResponse(double trueEnergy, 
                              const std::string& srcName) const {
   const std::string & diffuseComponent(diffuseSrcName(srcName));
   int indx(0);
   if (m_useEdisp) {
      indx = static_cast<int>((trueEnergy - m_trueEnergies[0])/m_estep);
      if (indx < 0 || indx >= static_cast<int>(m_trueEnergies.size())) {
         return 0;
      }
   }
   std::map<std::string, diffuse_response>::const_iterator it;
   if ((it = m_respDiffuseSrcs.find(diffuseComponent))
       != m_respDiffuseSrcs.end()) {
      if (m_useEdisp) {
         const diffuse_response & resp = it->second;
         double my_value = (trueEnergy - m_trueEnergies[indx])
            /(m_trueEnergies[indx+1] - m_trueEnergies[indx])
            *(resp[indx+1] - resp[indx]) + resp[indx];
         return my_value;
      } else {
// The response is just the single value in the diffuse_response vector.
         return it->second[0];
      }
   } else {
      std::string errorMessage 
         = "Event::diffuseResponse: \nDiffuse component " 
         + diffuseComponent 
         + " does not have an associated diffuse response.\n";
      throw Exception(errorMessage);
   }
   return 0;
}

const std::vector<double> & 
Event::diffuseResponse(const std::string& srcName) const {
   const std::string & diffuseComponent(diffuseSrcName(srcName));
   std::map<std::string, diffuse_response>::const_iterator it;
   if ((it = m_respDiffuseSrcs.find(diffuseComponent))
       == m_respDiffuseSrcs.end()) {
      std::string errorMessage 
         = "Event::diffuseResponse: \nDiffuse component " 
         + diffuseComponent
         + " does not have an associated diffuse response.\n";
      throw Exception(errorMessage);
   }
   return it->second;
}

void Event::computeResponseGQ(std::vector<DiffuseSource *> & srcList, 
                              const ResponseFunctions & respFuncs,
                              bool useDummyValue) {
   std::vector<DiffuseSource *> srcs;
   getNewDiffuseSrcs(srcList, srcs);
   if (srcs.size() == 0) {
      return;
   }
   double mumin(-1);
   double mumax(1);
   double err(1e-2);
   int ierr;

   EquinoxRotation eqRot(getDir().ra(), getDir().dec());
   for (size_t i(0); i < srcs.size(); i++) {
      const std::string& name(diffuseSrcName(srcs.at(i)->getName()));
      if (useDummyValue) {
         m_respDiffuseSrcs[name].push_back(0);
      } else {
         double respValue(0);
         mumin = -1;
         mumax = 1;
         double phimin(0);
         double phimax(2.*M_PI);
         try {
            srcs.at(i)->mapBaseObject()->getDiffRespLimits(getDir(), 
                                                           mumin, mumax,
                                                           phimin, phimax);
         } catch (MapBaseException &) {
            // do nothing
         }
         DiffRespIntegrand muIntegrand(*this, respFuncs, *srcs.at(i), eqRot,
                                       phimin, phimax);
         respValue = 
            st_facilities::GaussianQuadrature::dgaus8(muIntegrand, mumin,
                                                      mumax, err, ierr);
         m_respDiffuseSrcs[name].push_back(respValue);
      }
   }
}

void Event::computeResponse(std::vector<DiffuseSource *> &srcList, 
                            const ResponseFunctions & respFuncs,
                            double sr_radius, double sr_radius2) {
   std::vector<DiffuseSource *> srcs;
   getNewDiffuseSrcs(srcList, srcs);
   if (srcs.size() == 0) {
      return;
   }
   double ra0(m_appDir.ra());
   double dec0(m_appDir.dec());
   EquinoxRotation eqRot(ra0, dec0);
   if (!s_haveSourceRegionData) {
      prepareSrData(sr_radius, sr_radius2);
   }

   const std::vector<double> & muArray =
      LogNormalMuDist::instance()->muPoints(m_energy);

// Create a vector of srcDirs looping over the source region locations.
   std::vector<astro::SkyDir> srcDirs;
   for (unsigned int i = 0; i < muArray.size(); i++) {
      for (unsigned int j = 0; j < s_phi.size(); j++) {
         astro::SkyDir srcDir;
         getCelestialDir(s_phi[j], muArray[i], eqRot, srcDir);
         srcDirs.push_back(srcDir);
      }
   }
   std::vector<double>::iterator trueEnergy = m_trueEnergies.begin();
   for ( ; trueEnergy != m_trueEnergies.end(); ++trueEnergy) {

// Prepare the array of integrals over phi for passing to the 
// trapezoidal integrator for integration over mu.
      std::vector< std::vector<double> > mu_integrands;
      mu_integrands.resize(srcs.size());
      for (unsigned int i = 0; i < muArray.size(); i++) {

// Prepare phi-integrand arrays.
         std::vector< std::vector<double> > phi_integrands;
         phi_integrands.resize(srcs.size());
         for (unsigned int j = 0; j < s_phi.size(); j++) {
            int indx = i*s_phi.size() + j;
            astro::SkyDir & srcDir = srcDirs[indx];
            double inc = m_scDir.SkyDir::difference(srcDir)*180./M_PI;
            if (inc < 90.) {
               double totalResp = 
                  respFuncs.totalResponse(*trueEnergy, m_energy,
                                          m_scDir, m_scXDir, srcDir, m_appDir,
                                          m_type);
               for (unsigned int k = 0; k < srcs.size(); k++) {
                  double srcDist_val 
                     = srcs[k]->spatialDist(SkyDirArg(srcDir, *trueEnergy));
                  if (srcDist_val < 0) {
                     throw std::runtime_error("srcDist_val < 0");
                  }
                  phi_integrands[k].push_back(totalResp*srcDist_val);
               }
            } else {
               for (unsigned int k = 0; k < srcs.size(); k++)
                  phi_integrands[k].push_back(0);
            }
         }
         
// Perform the phi-integrals
         for (unsigned int k = 0; k < srcs.size(); k++) {
            TrapQuad phiQuad(s_phi, phi_integrands[k]);
            mu_integrands[k].push_back(phiQuad.integral());
         }
      }

// Perform the mu-integrals
      for (unsigned int k = 0; k < srcs.size(); k++) {
         TrapQuad muQuad(muArray, mu_integrands[k]);
	 const std::string & name(diffuseSrcName(srcs[k]->getName()));
         double respValue = muQuad.integral();
         if (respValue < 0) {
            throw std::runtime_error("Negative diffuse response value computed"
                                     " in Likelihood::Event::computeResponse");
         }
         m_respDiffuseSrcs[name].push_back(muQuad.integral());
      }
   } // loop over trueEnergy
// // Compute the Gaussian params to check validity of response.
//    for (unsigned int k = 0; k < srcs.size(); k++) {
//       std::string name = srcs[k]->getName();
//       double norm, mean, sigma;
//       computeGaussianParams(name, norm, mean, sigma);
//    }
}

void Event::computeGaussianParams(const std::string & name,
                                  double &norm, double &mean, 
                                  double &sigma) const {
   norm = 0;
   double eavg(0);
   double e2avg(0);
   const std::vector<double> & resp = diffuseResponse(name);
   for (unsigned int i = 0; i < m_trueEnergies.size() - 1; i++) {
      double deltaE = (m_trueEnergies[i+1] - m_trueEnergies[i]);
      norm += deltaE*(resp[i] + resp[i+1])/2.;
      eavg += deltaE*(m_trueEnergies[i]*resp[i] 
                      + m_trueEnergies[i+1]*resp[i+1])/2.;
      e2avg += deltaE*(m_trueEnergies[i]*m_trueEnergies[i]*resp[i] 
                       + m_trueEnergies[i+1]*m_trueEnergies[i+1]*resp[i+1])/2.;
   }
   if (norm <= 0) {    // should never happen...
      mean = m_energy;
      sigma = m_energy/10.;
   } else {
      mean = eavg/norm;
      sigma = sqrt(e2avg/norm - mean*mean);
   }
}

void Event::writeDiffuseResponses(const std::string & filename) {
   std::ofstream outfile(filename.c_str());
   std::map<std::string, diffuse_response>::iterator it
      = m_respDiffuseSrcs.begin();
   for ( ; it != m_respDiffuseSrcs.end(); ++it) {
      diffuse_response & resp = it->second;
      for (unsigned int ie = 0; ie < resp.size(); ie++) {
         outfile << m_trueEnergies[ie] << "  "
                 << resp[ie] << std::endl;
      }
   }
   outfile.close();
}

void Event::prepareSrData(double sr_radius, double sr_radius2) {
   fillMuArray(sr_radius, 100, s_mu);
   fillMuArray(sr_radius2, 200, s_mu_2);
   int nphi(50);
   double phistep = 2.*M_PI/(nphi - 1.);
   for (int i = 0; i < nphi; i++) {
      s_phi.push_back(phistep*i);
   }
   s_haveSourceRegionData = true;
}

void Event::fillMuArray(double sr_radius, int nmu, 
                        std::vector<double> & mu) const {
   (void)(sr_radius);
   (void)(nmu);
//    double mumin = cos(sr_radius*M_PI/180);
// // Sample more densely near theta = 0:
//    std::deque<double> my_mu;
//    double nscale = static_cast<double>((nmu-1)*(nmu-1));
//    for (int i = 0; i < nmu; i++) {
//       my_mu.push_front(1. - i*i/nscale*(1. - mumin));
//    }
//    mu.resize(my_mu.size());
//    std::copy(my_mu.begin(), my_mu.end(), mu.begin());

   double one_m_mu[] = {5.22906e-06, 1.13481e-05, 1.85455e-05, 2.68409e-05,
                        3.62683e-05, 4.68339e-05, 5.86167e-05, 7.16618e-05,
                        8.60305e-05, 1.01793e-04, 1.19020e-04, 1.37776e-04,
                        1.58140e-04, 1.80313e-04, 2.04256e-04, 2.30022e-04,
                        2.58057e-04, 2.88043e-04, 3.20447e-04, 3.55255e-04,
                        3.92521e-04, 4.32754e-04, 4.75544e-04, 5.21851e-04,
                        5.70961e-04, 6.24026e-04, 6.80406e-04, 7.40933e-04,
                        8.05721e-04, 8.74407e-04, 9.48961e-04, 1.02656e-03,
                        1.11241e-03, 1.20181e-03, 1.29860e-03, 1.40184e-03,
                        1.51031e-03, 1.62982e-03, 1.75387e-03, 1.88928e-03,
                        2.03329e-03, 2.18415e-03, 2.35182e-03, 2.52568e-03,
                        2.71457e-03, 2.91773e-03, 3.12821e-03, 3.36548e-03,
                        3.61242e-03, 3.87667e-03, 4.16746e-03, 4.46873e-03,
                        4.80398e-03, 5.16041e-03, 5.53458e-03, 5.95799e-03,
                        6.39695e-03, 6.87942e-03, 7.40369e-03, 7.94715e-03,
                        8.57575e-03, 9.22864e-03, 9.94570e-03, 1.07339e-02,
                        1.15550e-02, 1.25117e-02, 1.35061e-02, 1.46251e-02,
                        1.58408e-02, 1.71580e-02, 1.86531e-02, 2.02231e-02,
                        2.20740e-02, 2.40050e-02, 2.62824e-02, 2.86946e-02,
                        3.15273e-02, 3.45642e-02, 3.81605e-02, 4.20163e-02,
                        4.66868e-02, 5.17741e-02, 5.78414e-02, 6.47894e-02,
                        7.27099e-02, 8.24593e-02, 9.37943e-02, 1.07163e-01,
                        1.23816e-01, 1.44448e-01, 1.70216e-01, 2.03157e-01,
                        2.46342e-01, 3.04421e-01, 3.87505e-01, 5.10091e-01,
                        7.07397e-01, 1.05980e+00, 1.73826e+00};
   std::deque<double> my_mu;
   for (size_t i(0); i < 99; i++) {
      my_mu.push_front(1. - one_m_mu[i]);
   }
   mu.resize(my_mu.size());
   std::copy(my_mu.begin(), my_mu.end(), mu.begin());
}

void Event::getCelestialDir(double phi, double mu, 
                            EquinoxRotation & eqRot,
                            astro::SkyDir & dir) {
   double sp = sin(phi);
   double arg = mu/sqrt(1 - (1 - mu*mu)*sp*sp);
   double alpha;
   if (cos(phi) < 0) {
      alpha = 2*M_PI - my_acos(arg);
   } else {
      alpha = my_acos(arg);
   }
   double delta = asin(sqrt(1 - mu*mu)*sp);

// The direction in "Equinox rotated" coordinates
   astro::SkyDir indir(alpha*180/M_PI, delta*180/M_PI);

// Convert to the unrotated coordinate system (should probably use 
// Hep3Vector methods here instead).
   eqRot.do_rotation(indir, dir);
}

void Event::getNewDiffuseSrcs(const std::vector<DiffuseSource *> & srcList,
                              std::vector<DiffuseSource *> & srcs) const {
   for (std::vector<DiffuseSource *>::const_iterator it = srcList.begin();
        it != srcList.end(); ++it) {
      const std::string & name(diffuseSrcName((*it)->getName()));
      if (!m_respDiffuseSrcs.count(name)) {
         srcs.push_back(*it);
      }
   }
}

void Event::toLower(std::string & name) {
   for (std::string::iterator it = name.begin(); it != name.end(); ++it) {
      *it = std::tolower(*it);
   }
}

const std::string & Event::diffuseSrcName(const std::string & srcName) const {
   std::map<std::string, std::string>::iterator it = 
     m_diffSrcNames.find(srcName);
   if(it == m_diffSrcNames.end())
     {
       std::string name(m_respName + "__" + srcName);
       toLower(name);
       m_diffSrcNames[srcName] = name;
       return m_diffSrcNames[srcName];
     }
   return it->second;
}

void Event::updateModelSum(const Source & src, CachedResponse* cResp) {
   double contribution(src.fluxDensity(*this, cResp));
   const std::string & srcName(src.getName());
   std::map<std::string, double>::iterator fluxDensity;
   if ((fluxDensity = m_fluxDensities.find(srcName)) 
       != m_fluxDensities.end()) {
      m_modelSum -= fluxDensity->second;
   }
   m_modelSum += contribution;
   m_fluxDensities[srcName] = contribution;
}

void Event::deleteSource(const std::string & srcName) {
   m_modelSum -= m_fluxDensities[srcName];
   m_fluxDensities.erase(srcName);
}

void Event::resetModelSum() {
   m_modelSum = 0;
}

void Event::setDiffuseResponse(const std::string& srcName,
                               const std::vector<double> & gaussianParams) {
   const std::string & diffuseComponent(diffuseSrcName(srcName));
   static double sqrt2pi = sqrt(2.*M_PI);
   std::vector<double>::const_iterator energy = m_trueEnergies.begin();
   m_respDiffuseSrcs[diffuseComponent].clear();
   for ( ; energy != m_trueEnergies.end(); ++energy) {
      double value = gaussianParams[0]/sqrt2pi/gaussianParams[2]
         *exp(-(*energy - gaussianParams[1])*(*energy - gaussianParams[1])
              /gaussianParams[2]/gaussianParams[2]/2.);
      m_respDiffuseSrcs[diffuseComponent].push_back(value);
   }
}

} // namespace Likelihood
