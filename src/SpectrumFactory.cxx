/** @file SpectrumFactory.cxx
 * @brief Implementation for the SpectrumFactory class, which applies the
 * Prototype pattern to return clones of various spectral components 
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SpectrumFactory.cxx,v 1.3 2003/03/22 01:22:51 jchiang Exp $
 */

#include "Likelihood/SpectrumFactory.h"
#include <assert.h>

namespace Likelihood {

SpectrumFactory::~SpectrumFactory() {
   std::map<std::string, Function *>::iterator it = m_prototypes.begin();
   for (; it != m_prototypes.end(); it++)
      delete it->second;
}

void SpectrumFactory::addFunc(const std::string &name, Function* func, 
                              bool fromClone) {
   if (!m_prototypes.count(name)) {
      if (fromClone) {
         m_prototypes[name] = func->clone();
      } else {
         m_prototypes[name] = func;
      }
   } else {
      std::cerr << "SpectrumFactory: A Function named "
                << name << " already exists!" << std::endl;
      assert(false);
   }
}

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
