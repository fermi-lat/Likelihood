/**
 * @file ResponseFunctions.cxx
 * @brief Implementation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/ResponseFunctions.cxx,v 1.7 2004/06/01 04:27:00 jchiang Exp $
 */

#include "Likelihood/ScData.h"
#include "Likelihood/ResponseFunctions.h"

namespace Likelihood {

ResponseFunctions * ResponseFunctions::s_instance = 0;

std::map<unsigned int, latResponse::Irfs *> ResponseFunctions::s_respPtrs;

bool ResponseFunctions::s_useEdisp(false);
   
double ResponseFunctions::totalResponse(double time, 
                                        double energy, double appEnergy,
                                        const astro::SkyDir &srcDir,
                                        const astro::SkyDir &appDir, 
                                        int type) {
   ScData * scData = ScData::instance();
   astro::SkyDir zAxis = scData->zAxis(time);
   astro::SkyDir xAxis = scData->xAxis(time);
   
   double myResponse(0);
   std::map<unsigned int, latResponse::Irfs *>::iterator respIt 
      = instance()->begin();
   for ( ; respIt != instance()->end(); respIt++) {
      if (respIt->second->irfID() == type) {  
         latResponse::IPsf * psf = respIt->second->psf();
         latResponse::IAeff * aeff = respIt->second->aeff();
         double psf_val = psf->value(appDir, energy, srcDir, zAxis, xAxis);
         double aeff_val = aeff->value(energy, srcDir, zAxis, xAxis);
         if (s_useEdisp) {
            latResponse::IEdisp * edisp = respIt->second->edisp();
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

latResponse::Irfs * ResponseFunctions::respPtr(unsigned int i) {
   if (s_respPtrs.count(i)) {
      return s_respPtrs[i];
   } else {
      return 0;
   }
}

ResponseFunctions * ResponseFunctions::instance() {
   if (s_instance == 0) {
      s_instance = new ResponseFunctions();
   }
   return s_instance;
}

} // namespace Likelihood
