/** 
 * @file Source.cxx
 * @brief Source class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Source.cxx,v 1.2 2003/03/22 01:22:51 jchiang Exp $
 */

#include <vector>
#include <string>

#include "Likelihood/Source.h"

namespace Likelihood {

Source::Source(const Source &rhs) {
   m_name = rhs.m_name;

// relegate deep copy of m_functions to the subclasses

//  // make a deep copy
//     m_functions.clear();
//     FuncMap::const_iterator it = rhs.m_functions.begin();
//     for(; it != rhs.m_functions.end(); it++) {
//        Function *funcptr = (*it).second->clone();
//        m_functions[(*it).first] = funcptr;
//     }
}

} // namespace Likelihood
