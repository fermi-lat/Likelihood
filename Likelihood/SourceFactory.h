/** 
 * @file SourceFactory.h
 * @brief Declaration of SourceFactory class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceFactory.h,v 1.18 2004/11/17 00:02:12 jchiang Exp $
 */

#ifndef Likelihood_SourceFactory_h
#define Likelihood_SourceFactory_h

#include <map>
#include <string>

#include <xercesc/dom/DOM.hpp>

#include "xml/Dom.h"

#include "Likelihood/Source.h"
#include "Likelihood/Exception.h"

namespace optimizers {
   class FunctionFactory;
}

namespace Likelihood {

#ifndef SWIG
using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;
#endif //SWIG

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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceFactory.h,v 1.18 2004/11/17 00:02:12 jchiang Exp $
 *
 */
    
class SourceFactory {

public:

   SourceFactory(bool verbose=false);

   virtual ~SourceFactory();

   Source *create(const std::string &name) throw(Exception);

   /// Clients should almost always have fromClone = true; otherwise,
   /// the destructor will delete their Source, rather than a clone.
   void addSource(const std::string &name, Source* src, 
                  bool fromClone = true) throw(Exception);

   void replaceSource(Source* src, bool fromClone = true);

   void readXml(const std::string &xmlFile,
                optimizers::FunctionFactory&,
                bool requireExposure=true) throw(Exception);

   void fetchSrcNames(std::vector<std::string> &srcNames);

private:

   bool m_verbose;

   std::map<std::string, Source *> m_prototypes;

   bool m_requireExposure;

#ifndef SWIG
   Source *makePointSource(const DOMElement * spectrum,
                           const DOMElement * spatialModel,
                           optimizers::FunctionFactory & funcFactory);

   Source *makeDiffuseSource(const DOMElement * spectrum,
                             const DOMElement * spatialModel,
                             optimizers::FunctionFactory & funcFactory);

   void setSpectrum(Source *src, const DOMElement *spectrum,
                    optimizers::FunctionFactory & funcFactory);
#endif

};

} // namespace Likelihood

#endif // Likelihood_SourceFactory_h
