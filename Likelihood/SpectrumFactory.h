/** 
 * @file SpectrumFactory.h
 * @brief Declaration of SpectrumFactory class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SpectrumFactory.h,v 1.10 2003/08/06 20:52:03 jchiang Exp $
 */

#ifndef Likelihood_SpectrumFactory_h
#define Likelihood_SpectrumFactory_h

#include <string>
#include <map>

#include "optimizers/../src/PowerLaw.h"
#include "optimizers/../src/Gaussian.h"
#include "optimizers/../src/AbsEdge.h"

#include "Likelihood/Exception.h"

namespace optimizers {
   class Function;
}

namespace Likelihood {

/** 
 * @class SpectrumFactory
 *
 * @brief This class implements the Prototype pattern to return
 * clones of various spectral components for building Sources.
 *
 * The design of this class is based on the Factory template class
 * of Hippodraw.
 *
 * In future, this class *may* be made Singleton.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SpectrumFactory.h,v 1.10 2003/08/06 20:52:03 jchiang Exp $
 *
 */
    
class SpectrumFactory {
public:

   SpectrumFactory() {
      addFunc("PowerLaw", new optimizers::PowerLaw(), false);
      addFunc("Gaussian", new optimizers::Gaussian(), false);
      addFunc("AbsEdge", new optimizers::AbsEdge(), false);
   }

   virtual ~SpectrumFactory();

   //! Clients should almost always have fromClone = true, unless
   //! they explicitly pass a new Function pointer; otherwise,
   //! the destructor will delete their Function.
   void addFunc(const std::string &name, optimizers::Function* func, 
                bool fromClone = true) throw(optimizers::Exception);

   optimizers::Function *makeFunction(const std::string &name);

   void listFunctions();

private:

   std::map<std::string, optimizers::Function *> m_prototypes;

};

} // namespace Likelihood

#endif // Likelihood_SpectrumFactory_h
