/** 
 * @file Source.cxx
 * @brief Source class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Source.cxx,v 1.3 2003/06/11 17:08:04 jchiang Exp $
 */

#include "Likelihood/Source.h"

namespace Likelihood {

Source::Source(const Source &rhs) {
// Delegate deep copy of m_functions to the subclasses.
   m_name = rhs.m_name;
   m_useEdisp = rhs.m_useEdisp;
}

} // namespace Likelihood
