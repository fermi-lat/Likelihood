/**
 * @file Pixel.cxx
 * @brief Flyweight representation of pixels on the sky that provides
 * access to model counts and derivatives wrt model parameters.
 * @author J. Chiang
 * 
 * $Header$
 */

#include "Likelihood/ExposureCube.h"
#include "Likelihood/Pixel.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/SourceModel.h"

namespace Likelihood {

SourceModel * Pixel::s_srcModel = new SourceModel();

double Pixel::modelCounts(double emin, double emax) const {
   std::vector<std::string> srcNames;
   s_srcModel->getSrcNames(srcNames);

   double my_counts(0);

   for (unsigned int i = 0; i < srcNames.size(); i++) {
      Source * src = s_srcModel->getSource(srcNames[i]);
      for (int evtType = 0; evtType < 2; evtType++) {
         Aeff aeff1(src, m_dir, emin, evtType);
         double map_lower = ExposureCube::instance()->value(m_dir, aeff1);
         Aeff aeff2(src, m_dir, emax, evtType);
         double map_upper = ExposureCube::instance()->value(m_dir, aeff2);
         my_counts += (map_lower + map_upper)/2.*m_solidAngle*(emax - emin);
      }
   }
   return my_counts;
}

void Pixel::getFreeDerivs(double emin, double emax, 
                          std::vector<double> & derivs) const {
   derivs.resize(s_srcModel->getNumFreeParams(), 0);

   unsigned int iparam(0);

   const std::map<std::string, Source *> & srcMap = s_srcModel->sources();

   std::map<std::string, Source *>::const_iterator src = srcMap.begin();
   for ( ; src != srcMap.end(); ++src) {
      Source::FuncMap srcFuncs = src->second->getSrcFuncs();
      Source::FuncMap::const_iterator func = srcFuncs.begin();
      for ( ; func != srcFuncs.end(); ++func) {
         std::vector<std::string> names;
         func->second->getFreeParamNames(names);
         for (unsigned int i = 0; i < names.size(); i++) {
            for (int evtType = 0; evtType < 2; evtType++) {
               AeffDeriv aeff1(src->second, names[i], m_dir, emin, evtType);
               double map1 = ExposureCube::instance()->value(m_dir, aeff1);
               AeffDeriv aeff2(src->second, names[i], m_dir, emax, evtType);
               double map2 = ExposureCube::instance()->value(m_dir, aeff2);
               derivs.at(iparam) += (map1 + map2)/2.*m_solidAngle
                  *(emax - emin);
            }
            iparam++;
         }
      }
   }
}

Pixel::Aeff::Aeff(Source * src, const astro::SkyDir & appDir, 
                  double energy, int type)
   : m_src(src), m_appDir(appDir), m_energy(energy), m_type(type) {
   PointSource * ptsrc = dynamic_cast<PointSource *>(src);
   if (ptsrc == 0) {
      m_separation = 90.;
   } else {
      m_separation = ptsrc->getDir().difference(appDir)*180./M_PI;
   }
}

double Pixel::Aeff::operator()(double costheta) const {
   if (m_separation < 90.) {
      double inclination = acos(costheta)*180./M_PI;
      static double phi(0);
      return m_src->fluxDensity(inclination, phi, m_energy, m_separation,
                                m_type);
   }
   return 0;
}

Pixel::AeffDeriv::AeffDeriv(Source * src, const std::string & paramName,
                            const astro::SkyDir & appDir, 
                            double energy, int type)
   : m_src(src), m_paramName(paramName), m_appDir(appDir), 
     m_energy(energy), m_type(type) {
   PointSource * ptsrc = dynamic_cast<PointSource *>(src);
   if (ptsrc == 0) {
      m_separation = 90.;
   } else {
      m_separation = ptsrc->getDir().difference(appDir)*180./M_PI;
   }
}

double Pixel::AeffDeriv::operator()(double costheta) const {
   if (m_separation < 90.) {
      double inclination = acos(costheta)*180./M_PI;
      static double phi(0);
      return m_src->fluxDensityDeriv(inclination, phi, m_energy, m_separation,
                                     m_type, m_paramName);
   }
   return 0;
}

} // namespace Likelihood
