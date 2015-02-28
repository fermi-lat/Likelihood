/** 
 * @file SourceFactory.h
 * @brief Declaration of SourceFactory class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceFactory.h,v 1.27 2011/09/16 23:20:28 sfegan Exp $
 */

#ifndef Likelihood_SourceFactory_h
#define Likelihood_SourceFactory_h

#include <map>
#include <string>

#include <xercesc/dom/DOM.hpp>

#include "xmlBase/Dom.h"

#include "Likelihood/Exception.h"
#include "Likelihood/Observation.h"
#include "Likelihood/Source.h"

namespace st_stream {
   class StreamFormatter;
}

namespace optimizers {
   class Function;
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
 */
    
class SourceFactory {

public:

   SourceFactory(const Observation & observation, bool verbose=false);

   ~SourceFactory();

   Source * create(const std::string & name);
   Source * releaseSource(const std::string & name);

   /// Clients should almost always have fromClone = true; otherwise,
   /// the destructor will delete their Source, rather than a clone.
   void addSource(const std::string &name, Source * src, 
                  bool fromClone=true);

   void replaceSource(Source * src, bool fromClone=true);

   void readXml(const std::string & xmlFile,
                optimizers::FunctionFactory &,
                bool requireExposure=true,
                bool addPointSources=true,
                bool loadMaps=true);

   void fetchSrcNames(std::vector<std::string> & srcNames);

private:

   bool m_verbose;

   std::map<std::string, Source *> m_prototypes;

   bool m_requireExposure;

   const Observation & m_observation;

   st_stream::StreamFormatter * m_formatter;

#ifndef SWIG
   Source *makePointSource(const DOMElement * spectrum,
                           const DOMElement * spatialModel,
                           optimizers::FunctionFactory & funcFactory);

   Source *makeDiffuseSource(const DOMElement * spectrum,
                             const DOMElement * spatialModel,
                             optimizers::FunctionFactory & funcFactory,
                             bool loadMap=true);

   void setSpectrum(Source *src, const DOMElement *spectrum,
                    optimizers::FunctionFactory & funcFactory);
#endif

   void checkRoiDist(double ra, double dec) const;

   std::string m_currentSrcName;

   void addParamsToMultipleBPL(optimizers::Function * spec, 
                               const std::vector<DOMElement *> & params,
                               const Source * src) const;

   void addParamsToPiecewisePL(optimizers::Function * spec, 
                               const std::vector<DOMElement *> & params,
                               const Source * src) const;

};

} // namespace Likelihood

#endif // Likelihood_SourceFactory_h
