/**
 * @file XmlBuilder.cxx
 * @brief Concrete implementation that is shareable by subclasses.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/XmlBuilder.cxx,v 1.3 2004/11/11 20:23:11 jchiang Exp $
 */

#include "xmlBase/XmlParser.h"

#include "optimizers/Dom.h"

#include "Likelihood/XmlBuilder.h"

namespace Likelihood {

XmlBuilder::XmlBuilder() {
   m_parser = new xmlBase::XmlParser();
   m_doc = optimizers::Dom::createDocument();
}

XmlBuilder::~XmlBuilder() {
   m_doc->release();
   delete m_parser;
}

} //namespace Likelihood
