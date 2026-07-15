/**
 * @file XmlBuilder.h
 * @brief Base class for XML builders using RapidXML
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/XmlBuilder.h,v 1.7 2011/06/27 00:16:19 jchiang Exp $
 */

#ifndef Likelihood_XmlBuilder_h
#define Likelihood_XmlBuilder_h

#include "xmlBase/rapidxml.hpp"

namespace Likelihood {

class XmlBuilder {
public:
   XmlBuilder();
   virtual ~XmlBuilder();

protected:
   rapidxml::xml_document<>* m_doc;
   rapidxml::xml_node<>* m_srcLib;
};

} // namespace Likelihood

#endif // Likelihood_XmlBuilder_h
