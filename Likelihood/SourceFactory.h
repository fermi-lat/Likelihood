/** 
 * @file SourceFactory.h
 * @brief Declaration of SourceFactory class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceFactory.h,v 1.9 2003/06/11 17:08:02 jchiang Exp $
 */

#ifndef Likelihood_SourceFactory_h
#define Likelihood_SourceFactory_h

#include <string>
#include <map>
#include "Likelihood/Source.h"
#include "Likelihood/LikelihoodException.h"

namespace Likelihood {

/** 
 * @class SourceFactory
 *
 * @brief This class implements the Prototype pattern to return
 * clones of various gamma-ray Sources.
 *
 * The design of this class is based on the Factory template class
 * of Hippodraw.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceFactory.h,v 1.9 2003/06/11 17:08:02 jchiang Exp $
 *
 */
    
class SourceFactory {
public:

   SourceFactory();

   virtual ~SourceFactory();

   //! Clients should almost always have fromClone = true; 
   //! otherwise, the destructor will delete their Source, rather than 
   //! a clone.
   void addSource(const std::string &name, Source* src, 
                  bool fromClone = true) throw(LikelihoodException);

   void replaceSource(Source* src, bool fromClone = true);

   Source *makeSource(const std::string &name);

   void fetchSrcNames(std::vector<std::string> &srcNames);

private:

   std::map<std::string, Source *> m_prototypes;

};

} // namespace Likelihood

#endif // Likelihood_SourceFactory_h
