/** @file Source.cxx
 * @brief Source class implementation
 * @author J. Chiang
 *
 * $Header$
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
