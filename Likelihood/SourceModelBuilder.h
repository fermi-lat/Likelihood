/**
 * @file SourceModelBuilder.h
 * @brief Builder class for creating xml files of Source components.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceModelBuilder.h,v 1.1 2004/02/20 00:02:03 jchiang Exp $
 */

#ifndef Likelihood_SourceModelBuilder_h
#define Likelihood_SourceModelBuilder_h

#include "Likelihood/XmlBuilder.h"

namespace optimizers {
   class Function;
}

namespace Likelihood {

/**
 * @class SourceModelBuilder
 * @brief This class provides methods for writing the source
 * model information as xml.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceModelBuilder.h,v 1.1 2004/02/20 00:02:03 jchiang Exp $
 */

class SourceModelBuilder : public XmlBuilder {

public:

   SourceModelBuilder(const std::string &functionLibrary,
                      const std::string &srcLibTitle);
   
   virtual ~SourceModelBuilder();

   virtual void addSource(Source &src);

   virtual void write(std::string xmlFile);

private:

   typedef XERCES_CPP_NAMESPACE_QUALIFIER DOMElement DomElement;

   DomElement * m_srcLib;

   DomElement * likelihoodSource(Source & src);
   DomElement * spectralPart(Source & src);
   void addSpatialPart(DomElement * srcElt, Source & src);

};

} // namespace Likelihood

#endif // Likelihood_SourceModelBuilder_h
