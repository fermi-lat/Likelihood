/** 
 * @file Source.cxx
 * @brief Source class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Source.cxx,v 1.5 2005/03/01 01:06:55 jchiang Exp $
 */

#include "Likelihood/Source.h"

namespace Likelihood {

Source::Source() : m_name(""), m_srcType(""), m_useEdisp(false) {}

Source::Source(const Source &rhs) {
// Delegate the deep copy of m_functions to the subclasses.
   m_name = rhs.m_name;
   m_srcType = rhs.m_srcType;
   m_useEdisp = rhs.m_useEdisp;
}

} // namespace Likelihood
