/** @file SourceFactory.h
 * @brief Declaration of SourceFactory class
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef SourceFactory_h
#define SourceFactory_h

#include <string>
#include <map>

#include "Likelihood/Source.h"
#include "Likelihood/PointSource.h"

namespace Likelihood {
/** 
 * @class SourceFactory
 *
 * @brief This class implements the Prototype pattern to return
 * clones of various gamma-ray Sources.
 *
 * @author J. Chiang
 *    
 * $Header$
 *
 */
    
class SourceFactory {
public:

   SourceFactory();

   virtual ~SourceFactory() {}

   void addSource(const std::string &name, Source* src, 
                  bool fromClone = true);

   Source *makeSource(const std::string &name);

   void listSources();

private:

   std::map<std::string, Source *> m_prototypes;

};

} // namespace Likelihood

#endif // SourceFactory_h
