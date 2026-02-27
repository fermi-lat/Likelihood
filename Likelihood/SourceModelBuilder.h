/**
 * @file SourceModelBuilder.h
 * @brief Builder class for creating xml files of Source components.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceModelBuilder.h,v 1.2 2004/11/11 00:03:28 jchiang Exp $
 */

#ifndef Likelihood_SourceModelBuilder_h
#define Likelihood_SourceModelBuilder_h

#include "Likelihood/XmlBuilder.h"
#include "xmlbase/rapidxml.h"

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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceModelBuilder.h,v 1.2 2004/11/11 00:03:28 jchiang Exp $
 */


//#ifndef SWIG
//typedef XERCES_CPP_NAMESPACE_QUALIFIER DOMElement DomElement;
//typedef XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument DOMDocument;
//#endif //SWIG


class SourceModelBuilder : public XmlBuilder {

public:

   SourceModelBuilder(const std::string &functionLibrary,
                      const std::string &srcLibTitle);

   virtual ~SourceModelBuilder();

   virtual void addSourceModel(const SourceModel& srcModel);

   virtual void addSource(const Source &src);

   virtual void write(std::string xmlFile);

protected:


private:

   xmlbase::xml_node<> * m_srcLib;

   xmlbase::xml_node<> * likelihoodSource(const Source & src);
   xmlbase::xml_node<> * spectralPart(const Source & src);
   void addSpatialPart(xmlbase::xml_node<> * srcElt, const Source & src);
   void addComposite(xmlbase::xml_node<> * srcElt, const Source & src);
   void append_source(xmlbase::xml_node<> * parent, const Source &src);   
   void append_source_model(xmlbase::xml_node<> * parent, const SourceModel& srcModel);   

};

} // namespace Likelihood

#endif // Likelihood_SourceModelBuilder_h
