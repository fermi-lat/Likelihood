/**
 * @file XmlBuilder.cxx
 * @brief Concrete implementation that is shareable by subclasses.
 * @author J. Chiang
 *
 * $Header$
 */

#include "Likelihood/XmlBuilder.h"

namespace {
   DomDocument * createDocument() {
      DomDocument * doc = new DOM_Document();
      *doc = DOM_Document::createDocument();
      return doc;
   }
}

namespace Likelihood {

XmlBuilder::XmlBuilder() {
   m_parser = new xml::XmlParser();
   m_doc = ::createDocument();
}

XmlBuilder::~XmlBuilder() {
   delete m_doc;
   delete m_parser;
}

} //namespace Likelihood
