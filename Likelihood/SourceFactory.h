/** 
 * @file SourceFactory.h
 * @brief Declaration of SourceFactory class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceFactory.h,v 1.15 2004/02/20 22:51:06 jchiang Exp $
 */

#ifndef Likelihood_SourceFactory_h
#define Likelihood_SourceFactory_h

#include <string>
#include <map>

#include "xml/Dom.h"

#include "Likelihood/Source.h"
#include "Likelihood/Exception.h"

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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceFactory.h,v 1.15 2004/02/20 22:51:06 jchiang Exp $
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

   Source *makePointSource(const DomElement &spectrum,
                           const DomElement &spatialModel,
                           optimizers::FunctionFactory &funcFactory);

   Source *makeDiffuseSource(const DomElement &spectrum,
                             const DomElement &spatialModel,
                             optimizers::FunctionFactory &funcFactory);

   void setSpectrum(Source *src, const DomElement &spectrum,
                    optimizers::FunctionFactory &funcFactory);

};

} // namespace Likelihood

#endif // Likelihood_SourceFactory_h
