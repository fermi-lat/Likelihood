/** @file SpectrumFactory.h
 * @brief Declaration of SpectrumFactory class
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef SpectrumFactory_h
#define SpectrumFactory_h

#include <string>
#include <map>

#include "Likelihood/Function.h"
#include "PowerLaw.h"
#include "Gaussian.h"
#include "AbsEdge.h"

namespace Likelihood {
/** 
 * @class SpectrumFactory
 *
 * @brief This class implements the Prototype pattern to return
 * clones of various spectral components for building Sources.
 *
 * In future, this class *may* be made Singleton.
 *
 * @author J. Chiang
 *    
 * $Header$
 *
 */
    
class SpectrumFactory {
public:

   SpectrumFactory() {
      addFunc("PowerLaw", new PowerLaw());
      addFunc("Gaussian", new Gaussian());
      addFunc("AbsEdge", new AbsEdge());
   }

   virtual ~SpectrumFactory() {}

   void addFunc(const std::string &name, Function* func) {
      if (!m_prototypes.count(name)) m_prototypes[name] = func;
   }

   Function *makeFunction(const std::string &name);

   void listFunctions();

private:

   std::map<std::string, Function *> m_prototypes;

};

} // namespace Likelihood

#endif // SpectrumFactory_h
