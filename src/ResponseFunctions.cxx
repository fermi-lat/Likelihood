/**
 * @file ResponseFunctions.cxx
 * @brief Implementation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/ResponseFunctions.cxx,v 1.1 2003/10/22 04:30:33 jchiang Exp $
 */

#include "Likelihood/ResponseFunctions.h"

namespace Likelihood {

ResponseFunctions * ResponseFunctions::s_instance = 0;

std::map<unsigned int, latResponse::Irfs *> ResponseFunctions::s_respPtrs;

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
