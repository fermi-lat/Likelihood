/**
 * @file ResponseFunctions.cxx
 * @brief Implementation.
 * @author J. Chiang
 *
 * $Header$
 */

#include "Likelihood/ResponseFunctions.h"

namespace Likelihood {

ResponseFunctions * ResponseFunctions::s_instance = 0;

latResponse::Irfs * ResponseFunctions::respPtr(unsigned int i) {
   if (m_respPtrs.count(i)) {
      return m_respPtrs[i];
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
