/**
 * @file ResponseFunctions.cxx
 * @brief Implementation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/ResponseFunctions.cxx,v 1.22 2006/02/28 20:15:33 jchiang Exp $
 */

#include <sstream>
#include <stdexcept>

#include "astro/SkyDir.h"

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
                                        int type) const {
   bool foundResponse(false);
   double myResponse(0);
   std::map<unsigned int, irfInterface::Irfs *>::const_iterator respIt 
      = begin();
   for ( ; respIt != end(); respIt++) {
      if (respIt->second->irfID() == type) {
         foundResponse = true;
         irfInterface::IPsf * psf = respIt->second->psf();
         irfInterface::IAeff * aeff = respIt->second->aeff();
         double psf_val = psf->value(appDir, energy, srcDir, zAxis, xAxis);
         double aeff_val = aeff->value(energy, srcDir, zAxis, xAxis);
         if (m_useEdisp) {
            irfInterface::IEdisp * edisp = respIt->second->edisp();
            double edisp_val = edisp->value(appEnergy, energy, srcDir, 
                                            zAxis, xAxis);
            myResponse += psf_val*aeff_val*edisp_val;
         } else {
            myResponse += psf_val*aeff_val;
         }            
      }
   }
   if (!foundResponse) {
      throw std::runtime_error("Could not find appropriate response functions "
                               "for these event data.");
   }
   return myResponse;
}

double ResponseFunctions::totalResponse(double inclination, double phi,
                                        double energy, double appEnergy,
                                        double separation, int type) const {
   bool foundResponse(false);
   double myResponse(0);
   std::map<unsigned int, irfInterface::Irfs *>::const_iterator respIt 
      = begin();
   for ( ; respIt != end(); respIt++) {
      if (respIt->second->irfID() == type) {  
         foundResponse = true;
         irfInterface::IPsf * psf = respIt->second->psf();
         irfInterface::IAeff * aeff = respIt->second->aeff();
         double psf_val = psf->value(separation, energy, inclination, phi);
         double aeff_val = aeff->value(energy, inclination, phi);
         if (m_useEdisp) {
            irfInterface::IEdisp * edisp = respIt->second->edisp();
            double edisp_val = edisp->value(appEnergy, energy, 
                                            inclination, phi);
            myResponse += psf_val*aeff_val*edisp_val;
         } else {
            myResponse += psf_val*aeff_val;
         }            
      }
   }
   if (!foundResponse) {
      throw std::runtime_error("Could not find appropriate response functions "
                               "for these event data.");
   }
   return myResponse;
}   

irfInterface::Irfs * ResponseFunctions::respPtr(unsigned int i) const {
   if (m_respPtrs.count(i)) {
      return m_respPtrs.find(i)->second;
   } else {
      return 0;
   }
}

void ResponseFunctions::load(const std::string & respFuncs) {
//   irfLoader::Loader::go();
   irfLoader::Loader_go();
   irfInterface::IrfsFactory * myFactory 
      = irfInterface::IrfsFactory::instance();
      
   typedef std::map< std::string, std::vector<std::string> > respMap;
//   const respMap & responseIds = irfLoader::Loader::respIds();
   const respMap & responseIds = irfLoader::Loader_respIds();
   respMap::const_iterator it;
   if ( (it = responseIds.find(respFuncs)) != responseIds.end() ) {
      const std::vector<std::string> & resps = it->second;
      for (unsigned int i = 0; i < resps.size(); i++) {
         addRespPtr(i, myFactory->create(resps[i]));
      }
      setRespName(respFuncs);
   } else {
      std::ostringstream message;
      message << "Invalid response function choice: " << respFuncs << "\n"
              << "Valid choices are \n";
      for (respMap::const_iterator resp = responseIds.begin();
           resp != responseIds.end(); ++resp) {
         message << resp->first << "\n";
      }
      throw std::invalid_argument(message.str());
   }
}

} // namespace Likelihood
