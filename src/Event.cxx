/** 
 * @file Event.cxx
 * @brief Event class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/Event.cxx,v 1.84.2.1 2013/02/18 13:48:22 sfegan Exp $
 */

#include <cctype>

#include <algorithm>
#include <deque>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "st_facilities/GaussianQuadrature.h"

#include "Likelihood/DiffRespIntegrand.h"
#include "Likelihood/DiffRespIntegrand2.h"
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

#include "LogNormalMuDist.h"

namespace {
   std::string strip_front_back(std::string respName) {
// Strip off ::FRONT or ::BACK qualifiers.
      size_t pos(respName.find("::"));
      if (pos != std::string::npos) {
         respName = respName.substr(0, pos);
      }
      return respName;
   }
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

Event::Event() : m_appDir(0, 0), m_energy(0), m_arrTime(0), m_muZenith(0),
                 m_type(0), m_classLevel(0), m_scDir(0, 0), m_scXDir(0, 0),
                 m_useEdisp(false),
                 m_respName(""), m_modelSum(0), m_fluxDensities(), m_estep(0),
                 m_trueEnergies(),
                 m_efficiency(1), m_respDiffuseSrcs(), m_diffSrcNames() {}

Event::Event(double ra, double dec, double energy, double time, 
             const astro::SkyDir & scZAxis, const astro::SkyDir & scXAxis, 
             double muZenith, bool useEdisp, const std::string & respName,
             int type, double efficiency) 
   : m_appDir(astro::SkyDir(ra, dec)), m_energy(energy), m_arrTime(time),
     m_muZenith(muZenith), m_type(type), m_classLevel(0), 
     m_scDir(scZAxis), m_scXDir(scXAxis),
     m_useEdisp(useEdisp), m_respName(respName), m_modelSum(0),
     m_fluxDensities(), m_estep(0), m_trueEnergies(),
     m_efficiency(efficiency), m_respDiffuseSrcs(), m_diffSrcNames() {
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
   double minusone(-1);
   double one(1);
   double mumin(minusone);
   double mumax(one);

   EquinoxRotation eqRot(getDir().ra(), getDir().dec());
   for (size_t i(0); i < srcs.size(); i++) {
      const std::string& name(diffuseSrcName(srcs.at(i)->getName()));
      if (useDummyValue) {
         m_respDiffuseSrcs[name].push_back(0);
      } else {
         double respValue(0);
         mumin = minusone;
         mumax = one;
         double phimin(0);
         double phimax(2.*M_PI);
         try {
            srcs.at(i)->mapBaseObject()->getDiffRespLimits(getDir(), 
                                                           mumin, mumax,
                                                           phimin, phimax);
         } catch (MapBaseException &) {
            // do nothing
         }
	 if (srcs.at(i)->mapBasedIntegral() || 
             (::getenv("MAP_BASED_DIFFRSP") 
              && (mumin != minusone || mumax != one))) {
            respValue = srcs.at(i)->diffuseResponse(*this);
	 } else {
            if (::getenv("USE_OLD_DIFFRSP")) {
               /// Old integration scheme with the phi integral
               /// evaluated inside the theta integral.
               respValue = DiffRespIntegrand::
                  do2DIntegration(*this, respFuncs, *srcs.at(i), eqRot,
                                  mumin, mumax, phimin, phimax, 0.01, 0.1);
            } else {
               /// Steve's integration scheme with the theta integral
               /// evaluated inside the phi integral.  The produces
               /// much more accurate results.
               respValue = DiffRespIntegrand2::
                  do2DIntegration(*this, respFuncs, *srcs.at(i), eqRot,
                                  mumin, mumax, phimin, phimax, 0.001, 0.01);
            }
	 }
         m_respDiffuseSrcs[name].push_back(respValue);
      }
   }
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
   if (it == m_diffSrcNames.end()) {
      std::string name(::strip_front_back(m_respName) + "__" + srcName);
      toLower(name);
      m_diffSrcNames[srcName] = name;
      return m_diffSrcNames[srcName];
   }
   return it->second;
}

void Event::updateModelSum(const Source & src, CachedResponse* cResp) {
   double contribution(src.fluxDensity(*this, cResp)*efficiency());
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
   m_fluxDensities.clear();
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
