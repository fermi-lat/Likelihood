/**
 * @file XmlBuilder.cxx
 * @brief Concrete implementation that is shareable by subclasses.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/XmlBuilder.cxx,v 1.6 2006/03/16 06:20:06 jchiang Exp $
 */

#include "xmlBase/XmlParser.h"

#include "optimizers/Dom.h"

#include "Likelihood/XmlBuilder.h"

#include "Likelihood/XmlParser.h"

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
