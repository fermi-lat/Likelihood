/**
 * @file ResponseFunctions.cxx
 * @brief Implementation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/ResponseFunctions.cxx,v 1.19 2005/11/16 20:53:50 jchiang Exp $
 */

#include <stdexcept>

#include "astro/SkyDir.h"

#include "irfInterface/IrfsFactory.h"

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
   double myResponse(0);
   std::map<unsigned int, irfInterface::Irfs *>::const_iterator respIt 
      = begin();
   for ( ; respIt != end(); respIt++) {
      if (respIt->second->irfID() == type) {  
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
   return myResponse;
}

double ResponseFunctions::totalResponse(double inclination, double phi,
                                        double energy, double appEnergy,
                                        double separation, int type) const {
   double myResponse(0);
   std::map<unsigned int, irfInterface::Irfs *>::const_iterator respIt 
      = begin();
   for ( ; respIt != end(); respIt++) {
      if (respIt->second->irfID() == type) {  
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
   irfLoader::Loader::go();
   irfInterface::IrfsFactory * myFactory 
      = irfInterface::IrfsFactory::instance();
      
   typedef std::map< std::string, std::vector<std::string> > respMap;
   const respMap & responseIds = irfLoader::Loader::respIds();
   respMap::const_iterator it;
   if ( (it = responseIds.find(respFuncs)) != responseIds.end() ) {
      const std::vector<std::string> & resps = it->second;
      for (unsigned int i = 0; i < resps.size(); i++) {
         addRespPtr(i, myFactory->create(resps[i]));
      }
      setRespName(respFuncs);
   } else {
      throw std::invalid_argument("Invalid response function choice: "
                                  + respFuncs);
   }
}

} // namespace Likelihood
