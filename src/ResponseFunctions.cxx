/**
 * @file ResponseFunctions.cxx
 * @brief Implementation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/ResponseFunctions.cxx,v 1.39 2014/12/23 05:42:29 jchiang Exp $
 */

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include "astro/SkyDir.h"

#include "st_stream/StreamFormatter.h"

#include "irfInterface/IrfsFactory.h"

#undef ST_API_EXPORTS
#include "irfLoader/Loader.h"

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
                                        int type, double time) const {
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
      double psf_val(psf->value(appDir, energy, srcDir, zAxis, xAxis, time));
      double aeff_val(aeff->value(energy, srcDir, zAxis, xAxis, time));
      if (m_useEdisp) {
         throw std::runtime_error("Attempt to use energy dispersion "
                                  "handling in unbinned analysis.");
         // irfInterface::IEdisp * edisp(irfs->edisp());
         // double edisp_val(edisp->value(appEnergy, energy, srcDir, 
         //                               zAxis, xAxis, time));
         // myResponse = psf_val*aeff_val*edisp_val;
      } else {
         myResponse = psf_val*aeff_val;
      }
   }
   return myResponse;
}

double ResponseFunctions::totalResponse(double inclination, double phi,
                                        double energy, double appEnergy,
                                        double separation, int type,
                                        double time) const {
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
      double psf_val(psf->value(separation, energy, inclination, phi, time));
      double aeff_val(aeff->value(energy, inclination, phi, time));
      if (m_useEdisp) {
         throw std::runtime_error("Attempt to use energy dispersion "
                                  "handling in unbinned analysis.");
         // irfInterface::IEdisp * edisp(irfs->edisp());
         // double edisp_val(edisp->value(appEnergy, energy, inclination, phi,
         //                               time));
         // myResponse = psf_val*aeff_val*edisp_val;
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
                             const std::vector<unsigned int> & evtTypes) {
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
            addRespPtr(irfId, irfs);
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
   if (m_respPtrs.size() == 0) {
      throw std::runtime_error("No valid irfs loaded for the "
                               "specified event selection.");
   }
   st_stream::StreamFormatter formatter("ResponseFunctions", "load", 3);
   formatter.info() << "ResponseFunctions::load: IRF used: " 
                    << m_respName << "\n"
                    << "  event_types:";
   for (std::map<unsigned int, irfInterface::Irfs *>::const_iterator 
           it(m_respPtrs.begin()); it != m_respPtrs.end(); ++ it) {
      formatter.info() << "  " << it->first;
   }
   formatter.info() << std::endl;
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
                                int type, double time) const {
   irfInterface::Irfs * irfs(const_cast<irfInterface::Irfs *>(respPtr(type)));
   if (!irfs) {
      std::ostringstream message;
      message << "Could not find appropriate response functions "
              << "for these event data."
              << "Event class requested: " << type << std::endl;
      throw std::runtime_error(message.str());
   }
   irfInterface::IEdisp * edisp(irfs->edisp());
   return edisp->value(emeas, etrue, appDir, zAxis, xAxis, time);
}

double ResponseFunctions::edisp(double emeas, double etrue, 
                                double theta, double phi, 
                                int type, double time) const {
   irfInterface::Irfs * irfs(const_cast<irfInterface::Irfs *>(respPtr(type)));
   if (!irfs) {
      std::ostringstream message;
      message << "Could not find appropriate response functions "
              << "for these event data."
              << "Event class requested: " << type << std::endl;
      throw std::runtime_error(message.str());
   }
   irfInterface::IEdisp * edisp(irfs->edisp());
   return edisp->value(emeas, etrue, theta, phi, time);
}

double ResponseFunctions::aeff(double etrue, 
                               const astro::SkyDir & appDir,
                               const astro::SkyDir & zAxis,
                               const astro::SkyDir & xAxis,
                               int type, double time) const {
   irfInterface::Irfs * irfs(const_cast<irfInterface::Irfs *>(respPtr(type)));
   if (!irfs) {
      std::ostringstream message;
      message << "Could not find appropriate response functions "
              << "for these event data."
              << "Event class requested: " << type << std::endl;
      throw std::runtime_error(message.str());
   }
   irfInterface::IAeff * aeff(irfs->aeff());
   return aeff->value(etrue, appDir, zAxis, xAxis, time);
}

double ResponseFunctions::aeff(double etrue, double theta, double phi,
                               int type, double time) const {
   irfInterface::Irfs * irfs(const_cast<irfInterface::Irfs *>(respPtr(type)));
   if (!irfs) {
      std::ostringstream message;
      message << "Could not find appropriate response functions "
              << "for these event data."
              << "Event class requested: " << type << std::endl;
      throw std::runtime_error(message.str());
   }
   irfInterface::IAeff * aeff(irfs->aeff());
   return aeff->value(etrue, theta, phi, time);
}

} // namespace Likelihood
