/**
 * @file XmlBuilder.cpp
 * @brief Concrete implementation that is shareable by subclasses.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/XmlBuilder.cxx,v 1.7 2011/06/27 00:16:19 jchiang Exp $
 */

#include "xmlBase/rapidxml.hpp"

#include "Likelihood/XmlBuilder.h"
#include "Likelihood/XmlParser.h"

namespace Likelihood {

XmlBuilder::XmlBuilder() {
   // Create a new XML document
   m_doc = new rapidxml::xml_document<>();
}

XmlBuilder::~XmlBuilder() {
   // RapidXML documents clean up their memory pools automatically when destroyed
   delete m_doc;
}

} // namespace Likelihood
