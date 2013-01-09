/**
 * @file Pixel.cxx
 * @brief Fly-weight representation of pixels on the sky that provides
 * access to model counts and derivatives wrt model parameters.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/Pixel.cxx,v 1.10 2012/03/28 22:00:43 jchiang Exp $
 */

#include "Likelihood/Pixel.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/SourceModel.h"

namespace Likelihood {

double Pixel::modelCounts(double emin, double emax, 
                          SourceModel & srcModel) const {
   std::vector<std::string> srcNames;
   srcModel.getSrcNames(srcNames);

   double my_counts(0);

   double epoch((srcModel.observation().expCube().tstart() 
                 + srcModel.observation().expCube().tstop())/2.);

   for (unsigned int i = 0; i < srcNames.size(); i++) {
      Source * src = srcModel.getSource(srcNames[i]);
      std::map<unsigned int, irfInterface::Irfs *>::const_iterator respIt 
         = srcModel.observation().respFuncs().begin();
      for ( ; respIt != srcModel.observation().respFuncs().end(); ++respIt) {
         int evtType = respIt->second->irfID();
         Aeff aeff1(src, m_dir, emin, evtType, epoch);
         double map_lower = 
            srcModel.observation().expCube().value(m_dir, aeff1, emin);
         Aeff aeff2(src, m_dir, emax, evtType, epoch);
         double map_upper = 
            srcModel.observation().expCube().value(m_dir, aeff2, emax);
//         my_counts += (map_lower + map_upper)/2.*m_solidAngle*(emax - emin);
         my_counts += ( (map_lower*emin + map_upper*emax)/2.*m_solidAngle
                        *std::log(emax/emin) );
      }
   }
   return my_counts;
}

void Pixel::getFreeDerivs(double emin, double emax, SourceModel & srcModel, 
                          std::vector<double> & derivs) const {
   derivs.resize(srcModel.getNumFreeParams(), 0);

   unsigned int iparam(0);

   double epoch((srcModel.observation().expCube().tstart() 
                 + srcModel.observation().expCube().tstop())/2.);

   const std::map<std::string, Source *> & srcMap = srcModel.sources();

   std::map<std::string, Source *>::const_iterator src = srcMap.begin();
   for ( ; src != srcMap.end(); ++src) {
      std::vector<std::string> names;
      src->second->spectrum().getFreeParamNames(names);
      for (size_t i(0); i < names.size(); i++) {
         std::map<unsigned int, irfInterface::Irfs *>::const_iterator 
            respIt = srcModel.observation().respFuncs().begin();
         for ( ; respIt != srcModel.observation().respFuncs().end();
               ++respIt) {
            int evtType = respIt->second->irfID();
            AeffDeriv aeff1(src->second, names[i], m_dir, emin,
                            evtType, epoch);
            double map1 = 
               srcModel.observation().expCube().value(m_dir, aeff1, emin);
            AeffDeriv aeff2(src->second, names[i], m_dir, emax,
                            evtType, epoch);
            double map2 = 
               srcModel.observation().expCube().value(m_dir, aeff2, emax);
            derivs.at(iparam) += (map1 + map2)/2.*m_solidAngle
               *(emax - emin);
         }
         iparam++;
      }
   }
}

Pixel::Aeff::Aeff(Source * src, const astro::SkyDir & appDir, 
                  double energy, int type, double time)
   : ExposureCube::AeffBase(), m_src(src), m_appDir(appDir), 
     m_energy(energy), m_type(type), m_time(time) {
   PointSource * ptsrc = dynamic_cast<PointSource *>(src);
   if (ptsrc == 0) {
      m_separation = 90.;
   } else {
      m_separation = ptsrc->getDir().difference(appDir)*180./M_PI;
   }
}

double Pixel::Aeff::value(double costheta, double phi) const {
   if (m_separation < 90.) {
      double inclination = acos(costheta)*180./M_PI;
      return m_src->fluxDensity(inclination, phi, m_energy, m_appDir,
                                m_type, m_time);
   }
   return 0;
}

Pixel::AeffDeriv::AeffDeriv(Source * src, const std::string & paramName,
                            const astro::SkyDir & appDir, 
                            double energy, int type, double time)
   : ExposureCube::AeffBase(), m_src(src), m_paramName(paramName), 
     m_appDir(appDir), m_energy(energy), m_type(type), m_time(time) {
   PointSource * ptsrc = dynamic_cast<PointSource *>(src);
   if (ptsrc == 0) {
      m_separation = 90.;
   } else {
      m_separation = ptsrc->getDir().difference(appDir)*180./M_PI;
   }
}

double Pixel::AeffDeriv::value(double costheta, double phi) const {
   if (m_separation < 90.) {
      double inclination = acos(costheta)*180./M_PI;
      return m_src->fluxDensityDeriv(inclination, phi, m_energy, m_appDir,
                                     m_type, m_time, m_paramName);
   }
   return 0;
}

} // namespace Likelihood
