/**
 * @file ResponseFunctions.cxx
 * @brief Implementation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/ResponseFunctions.cxx,v 1.33 2011/10/17 05:54:07 jchiang Exp $
 */

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include "astro/SkyDir.h"

#include "irfInterface/IrfsFactory.h"

#undef ST_API_EXPORTS
#include "irfLoader/Loader.h"

#include "latResponse/IrfLoader.h"

#include "Likelihood/ResponseFunctions.h"

namespace Likelihood {
   
ResponseFunctions::~ResponseFunctions() {
   std::map<unsigned int, irfInterface::Irfs *>::iterator 
      it(m_respPtrs.begin());
   for ( ; it != m_respPtrs.end(); ++it) {
      deleteRespPtr(it->first);
   }
}

double ResponseFunctions::totalResponse(double energy, double appEnergy,
                                        const astro::SkyDir & zAxis,
                                        const astro::SkyDir & xAxis,
                                        const astro::SkyDir & srcDir,
                                        const astro::SkyDir & appDir, 
                                        int type) const {
   double myResponse;
   irfInterface::Irfs * irfs(const_cast<irfInterface::Irfs *>(respPtr(type)));
   if (!irfs) {
      std::ostringstream message;
      message << "Could not find appropriate response functions "
              << "for these event data."
              << "Event class requested: " << type << std::endl;
      throw std::runtime_error(message.str());
   } else {
      irfInterface::IPsf * psf(irfs->psf());
      irfInterface::IAeff * aeff(irfs->aeff());
      double psf_val(psf->value(appDir, energy, srcDir, zAxis, xAxis));
      double aeff_val(aeff->value(energy, srcDir, zAxis, xAxis));
      if (m_useEdisp) {
         irfInterface::IEdisp * edisp(irfs->edisp());
         double edisp_val(edisp->value(appEnergy, energy, srcDir, 
                                       zAxis, xAxis));
         myResponse = psf_val*aeff_val*edisp_val;
      } else {
         myResponse = psf_val*aeff_val;
      }
   }
   return myResponse;
}

double ResponseFunctions::totalResponse(double inclination, double phi,
                                        double energy, double appEnergy,
                                        double separation, int type) const {
   double myResponse;
   irfInterface::Irfs * irfs(const_cast<irfInterface::Irfs *>(respPtr(type)));
   if (!irfs) {
      std::ostringstream message;
      message << "Could not find appropriate response functions "
              << "for these event data." << std::endl
              << "Event class requested: " << type << std::endl;
      throw std::runtime_error(message.str());
   } else {
      irfInterface::IPsf * psf(irfs->psf());
      irfInterface::IAeff * aeff(irfs->aeff());
      double psf_val(psf->value(separation, energy, inclination, phi));
      double aeff_val(aeff->value(energy, inclination, phi));
      if (m_useEdisp) {
         irfInterface::IEdisp * edisp(irfs->edisp());
         double edisp_val(edisp->value(appEnergy, energy, inclination, phi));
         myResponse = psf_val*aeff_val*edisp_val;
      } else {
         myResponse = psf_val*aeff_val;
      }
   }
   return myResponse;
}   

irfInterface::Irfs * ResponseFunctions::respPtr(unsigned int i) const {
   typedef std::map<unsigned int, irfInterface::Irfs *> IrfMap_t;
   IrfMap_t::const_iterator irfs(m_respPtrs.find(i));
   if (irfs != m_respPtrs.end()) {
      return irfs->second;
   } else {
      return 0;
   }
}

void ResponseFunctions::load(const std::string & respFuncs,
                             const std::string & respBase,
                             const std::vector<size_t> & evtTypes) {
   // Turn of energy dispersion interpolation in
   // irfs/latResponse/Edisp2 for greater efficiency.
   latResponse::IrfLoader::set_edisp_interpolation(true);

   irfLoader::Loader_go();

   irfInterface::IrfsFactory * myFactory(irfInterface::IrfsFactory::instance());
      
   typedef std::map< std::string, std::vector<std::string> > respMap_t;
   const respMap_t & responseIds = irfLoader::Loader_respIds();
   respMap_t::const_iterator it;
   if ( (it = responseIds.find(respFuncs)) != responseIds.end() ) {
      const std::vector<std::string> & resps = it->second;
      for (unsigned int i = 0; i < resps.size(); i++) {
         irfInterface::Irfs * irfs(myFactory->create(resps[i]));
         size_t irfId(irfs->irfID());
         if (evtTypes.empty() ||
             std::count(evtTypes.begin(), evtTypes.end(), irfId) > 0) {
            addRespPtr(irfs->irfID(), irfs);
         } else {
            delete irfs;
         }
      }
      if (respBase == "") {
         setRespName(respFuncs);
      } else {
         setRespName(respBase);
      }
   } else {
      std::ostringstream message;
      message << "Invalid response function choice: " << respFuncs << "\n"
              << "Valid choices are \n";
      for (respMap_t::const_iterator resp = responseIds.begin();
           resp != responseIds.end(); ++resp) {
         message << resp->first << "\n";
      }
      throw std::invalid_argument(message.str());
   }
}

irfInterface::IEdisp & ResponseFunctions::edisp(int type) const {
   irfInterface::Irfs * irfs(const_cast<irfInterface::Irfs *>(respPtr(type)));
   if (!irfs) {
      std::ostringstream message;
      message << "Could not find appropriate response functions "
              << "for these event data."
              << "Event class requested: " << type << std::endl;
      throw std::runtime_error(message.str());
   }
   return *irfs->edisp();
}

double ResponseFunctions::edisp(double emeas, double etrue, 
                                const astro::SkyDir & appDir,
                                const astro::SkyDir & zAxis,
                                const astro::SkyDir & xAxis,
                                int type) const {
   irfInterface::Irfs * irfs(const_cast<irfInterface::Irfs *>(respPtr(type)));
   if (!irfs) {
      std::ostringstream message;
      message << "Could not find appropriate response functions "
              << "for these event data."
              << "Event class requested: " << type << std::endl;
      throw std::runtime_error(message.str());
   }
   irfInterface::IEdisp * edisp(irfs->edisp());
   return edisp->value(emeas, etrue, appDir, zAxis, xAxis);
}

double ResponseFunctions::edisp(double emeas, double etrue, 
                                double theta, double phi, int type) const {
   irfInterface::Irfs * irfs(const_cast<irfInterface::Irfs *>(respPtr(type)));
   if (!irfs) {
      std::ostringstream message;
      message << "Could not find appropriate response functions "
              << "for these event data."
              << "Event class requested: " << type << std::endl;
      throw std::runtime_error(message.str());
   }
   irfInterface::IEdisp * edisp(irfs->edisp());
   return edisp->value(emeas, etrue, theta, phi);
}

double ResponseFunctions::aeff(double etrue, 
                               const astro::SkyDir & appDir,
                               const astro::SkyDir & zAxis,
                               const astro::SkyDir & xAxis,
                               int type) const {
   irfInterface::Irfs * irfs(const_cast<irfInterface::Irfs *>(respPtr(type)));
   if (!irfs) {
      std::ostringstream message;
      message << "Could not find appropriate response functions "
              << "for these event data."
              << "Event class requested: " << type << std::endl;
      throw std::runtime_error(message.str());
   }
   irfInterface::IAeff * aeff(irfs->aeff());
   return aeff->value(etrue, appDir, zAxis, xAxis);
}

double ResponseFunctions::aeff(double etrue, double theta, double phi,
                               int type) const {
   irfInterface::Irfs * irfs(const_cast<irfInterface::Irfs *>(respPtr(type)));
   if (!irfs) {
      std::ostringstream message;
      message << "Could not find appropriate response functions "
              << "for these event data."
              << "Event class requested: " << type << std::endl;
      throw std::runtime_error(message.str());
   }
   irfInterface::IAeff * aeff(irfs->aeff());
   return aeff->value(etrue, theta, phi);
}

} // namespace Likelihood
