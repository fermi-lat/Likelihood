/**
 * @file XmlBuilder.cxx
 * @brief Concrete implementation that is shareable by subclasses.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/XmlBuilder.cxx,v 1.7 2011/06/27 00:16:19 jchiang Exp $
 */

#include "xmlBase/XmlParser.h"

#include "optimizers/Dom.h"

#include "Likelihood/XmlBuilder.h"

#include "Likelihood/XmlParser.h"

namespace Likelihood {

XmlBuilder::XmlBuilder()  {
//   m_parser = XmlParser::instance();
   m_parser = XmlParser_instance();
   m_doc = optimizers::Dom::createDocument();
}


XmlBuilder::~XmlBuilder() {
   m_doc->release();
}

} //namespace Likelihood
