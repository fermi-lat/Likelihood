/** @file SpectrumFactory.cxx
 * @brief Implementation for the SpectrumFactory class, which applies the
 * Prototype pattern to return clones of various spectral components 
 * @author J. Chiang
 *
 * $Header$
 */

#include "Likelihood/SpectrumFactory.h"

namespace Likelihood {

Function *SpectrumFactory::makeFunction(const std::string &name) {
   assert(m_prototypes.count(name));
   return m_prototypes[name]->clone();
}

void SpectrumFactory::listFunctions() {
   std::cout << "SpectrumFactory Functions: " << std::endl;
   std::map<std::string, Function *>::const_iterator 
      it = m_prototypes.begin();
   for (; it != m_prototypes.end(); it++)
      std::cout << it->first << std::endl;
}

} // namespace Likelihood
