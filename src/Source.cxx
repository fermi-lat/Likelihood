/** 
 * @file Source.cxx
 * @brief Source class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Source.cxx,v 1.4 2005/02/28 18:38:46 jchiang Exp $
 */

#include "Likelihood/Source.h"

namespace Likelihood {

Source::Source() {}

Source::Source(const Source &rhs) {
// Delegate deep copy of m_functions to the subclasses.
   m_name = rhs.m_name;
   m_useEdisp = rhs.m_useEdisp;
}

} // namespace Likelihood
