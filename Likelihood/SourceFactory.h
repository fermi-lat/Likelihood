/** 
 * @file SourceFactory.h
 * @brief Declaration of SourceFactory class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceFactory.h,v 1.12 2003/08/06 20:52:03 jchiang Exp $
 */

#ifndef Likelihood_SourceFactory_h
#define Likelihood_SourceFactory_h

#include <string>
#include <map>
#include "Likelihood/Source.h"
#include "Likelihood/Exception.h"

class DOM_Element;

namespace optimizers {

class FunctionFactory;

}

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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceFactory.h,v 1.12 2003/08/06 20:52:03 jchiang Exp $
 *
 */
    
class SourceFactory {

public:

   SourceFactory();

   virtual ~SourceFactory();

   Source *create(const std::string &name) throw(Exception);

   /// Clients should almost always have fromClone = true; otherwise,
   /// the destructor will delete their Source, rather than a clone.
   void addSource(const std::string &name, Source* src, 
                  bool fromClone = true) throw(Exception);

   void replaceSource(Source* src, bool fromClone = true);

   void readXml(const std::string &xmlFile,
                optimizers::FunctionFactory&) throw(Exception);

   void fetchSrcNames(std::vector<std::string> &srcNames);

private:

   std::map<std::string, Source *> m_prototypes;

   Source *makePointSource(const DOM_Element &spectrum,
                           const DOM_Element &spatialModel,
                           optimizers::FunctionFactory &funcFactory);

   Source *makeDiffuseSource(const DOM_Element &spectrum,
                             const DOM_Element &spatialModel,
                             optimizers::FunctionFactory &funcFactory);

   void setSpectrum(Source *src, const DOM_Element &spectrum,
                    optimizers::FunctionFactory &funcFactory);

};

} // namespace Likelihood

#endif // Likelihood_SourceFactory_h
