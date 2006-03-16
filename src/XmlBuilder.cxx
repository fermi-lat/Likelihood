/**
 * @file XmlBuilder.cxx
 * @brief Concrete implementation that is shareable by subclasses.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/XmlBuilder.cxx,v 1.5 2005/11/16 20:00:32 jchiang Exp $
 */

#include "xmlBase/XmlParser.h"

#include "optimizers/Dom.h"

#include "Likelihood/XmlBuilder.h"

#include "XmlParser.h"

namespace Likelihood {

XmlBuilder::XmlBuilder() {
//   m_parser = XmlParser::instance();
   m_parser = XmlParser_instance();
   m_doc = optimizers::Dom::createDocument();
}

XmlBuilder::~XmlBuilder() {
   m_doc->release();
}

} //namespace Likelihood
