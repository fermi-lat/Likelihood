/** 
 * @file SpectrumFactory.h
 * @brief Declaration of SpectrumFactory class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SpectrumFactory.h,v 1.4 2003/05/29 20:10:30 jchiang Exp $
 */

#ifndef SpectrumFactory_h
#define SpectrumFactory_h

#include <string>
#include <map>

//#include "Likelihood/Function.h"
#include "PowerLaw.h"
#include "Gaussian.h"
#include "AbsEdge.h"
#include "Likelihood/LikelihoodException.h"

namespace Likelihood {

class Function;

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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SpectrumFactory.h,v 1.4 2003/05/29 20:10:30 jchiang Exp $
 *
 */
    
class SpectrumFactory {
public:

   SpectrumFactory() {
      addFunc("PowerLaw", new PowerLaw(), false);
      addFunc("Gaussian", new Gaussian(), false);
      addFunc("AbsEdge", new AbsEdge(), false);
   }

   virtual ~SpectrumFactory();

   //! Clients should almost always have fromClone = true, unless
   //! they explicitly pass a new Function pointer; otherwise,
   //! the destructor will delete their Function.
   void addFunc(const std::string &name, Function* func, 
                bool fromClone = true) throw(LikelihoodException);

   Function *makeFunction(const std::string &name);

   void listFunctions();

private:

   std::map<std::string, Function *> m_prototypes;

};

} // namespace Likelihood

#endif // SpectrumFactory_h
