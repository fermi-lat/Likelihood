/**
 * @file ResponseFunctions.cxx
 * @brief Implementation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/ResponseFunctions.cxx,v 1.12 2004/12/01 16:46:27 jchiang Exp $
 */

#include "Likelihood/ScData.h"
#include "Likelihood/ResponseFunctions.h"

namespace Likelihood {

ResponseFunctions * ResponseFunctions::s_instance = 0;

std::map<unsigned int, irfInterface::Irfs *> ResponseFunctions::s_respPtrs;

bool ResponseFunctions::s_useEdisp(false);

std::string ResponseFunctions::s_respName("");
   
double ResponseFunctions::totalResponse(double time, 
                                        double energy, double appEnergy,
                                        const astro::SkyDir &srcDir,
                                        const astro::SkyDir &appDir, 
                                        int type) {
   ScData * scData = ScData::instance();
   astro::SkyDir zAxis = scData->zAxis(time);
   astro::SkyDir xAxis = scData->xAxis(time);
   
   double myResponse(0);
   std::map<unsigned int, irfInterface::Irfs *>::iterator respIt 
      = instance()->begin();
   for ( ; respIt != instance()->end(); respIt++) {
      if (respIt->second->irfID() == type) {  
         irfInterface::IPsf * psf = respIt->second->psf();
         irfInterface::IAeff * aeff = respIt->second->aeff();
         double psf_val = psf->value(appDir, energy, srcDir, zAxis, xAxis);
         double aeff_val = aeff->value(energy, srcDir, zAxis, xAxis);
         if (s_useEdisp) {
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
   
double ResponseFunctions::totalResponse(double energy, double appEnergy,
                                        const astro::SkyDir & zAxis,
                                        const astro::SkyDir & xAxis,
                                        const astro::SkyDir & srcDir,
                                        const astro::SkyDir & appDir, 
                                        int type) {
   double myResponse(0);
   std::map<unsigned int, irfInterface::Irfs *>::iterator respIt 
      = instance()->begin();
   for ( ; respIt != instance()->end(); respIt++) {
      if (respIt->second->irfID() == type) {  
         irfInterface::IPsf * psf = respIt->second->psf();
         irfInterface::IAeff * aeff = respIt->second->aeff();
         double psf_val = psf->value(appDir, energy, srcDir, zAxis, xAxis);
         double aeff_val = aeff->value(energy, srcDir, zAxis, xAxis);
         if (s_useEdisp) {
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
                                        double separation, int type) {
   double myResponse(0);
   std::map<unsigned int, irfInterface::Irfs *>::iterator respIt 
      = instance()->begin();
   for ( ; respIt != instance()->end(); respIt++) {
      if (respIt->second->irfID() == type) {  
         irfInterface::IPsf * psf = respIt->second->psf();
         irfInterface::IAeff * aeff = respIt->second->aeff();
         double psf_val = psf->value(separation, energy, inclination, phi);
         double aeff_val = aeff->value(energy, inclination, phi);
         if (s_useEdisp) {
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

irfInterface::Irfs * ResponseFunctions::respPtr(unsigned int i) {
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
