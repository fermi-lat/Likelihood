/**
 * @file XmlBuilder.cxx
 * @brief Concrete implementation that is shareable by subclasses.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/XmlBuilder.cxx,v 1.1 2004/02/20 00:02:07 jchiang Exp $
 */

#include "xml/XmlParser.h"

#include "optimizers/Dom.h"

#include "Likelihood/XmlBuilder.h"

namespace Likelihood {

XmlBuilder::XmlBuilder() {
   m_parser = new xml::XmlParser();
   m_doc = optimizers::Dom::createDocument();
}

XmlBuilder::~XmlBuilder() {
   delete m_doc;
   delete m_parser;
}

} //namespace Likelihood
