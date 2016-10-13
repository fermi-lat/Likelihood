/**
 * @file XmlBuilder.h
 * @brief Interface to builder classes for creating xml documents for
 * the flux and Likelihood packages.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/XmlBuilder.h,v 1.3 2005/01/03 23:01:16 jchiang Exp $
 */

#ifndef Likelihood_XmlBuilder_h
#define Likelihood_XmlBuilder_h

#include <xercesc/util/XercesDefs.hpp>

#include <string>

namespace xmlBase {
   class XmlParser;
}

namespace Likelihood {

   class Source;
   class SourceModel;

/**
 * @class XmlBuilder
 *
 * @brief Base class and interface definition for xml document
 * builders.  The Builder pattern traditionally provides empty default
 * implementations for all methods rather than making them pure
 * virtual [GOF].
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/XmlBuilder.h,v 1.3 2005/01/03 23:01:16 jchiang Exp $
 */

#ifndef SWIG
typedef XERCES_CPP_NAMESPACE_QUALIFIER DOMElement DOMElement;
typedef XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument DOMDocument;
#endif // SWIG

class XmlBuilder {

public:

   virtual ~XmlBuilder();

   virtual void addSourceModel(const SourceModel& srcModel) {}

   virtual void addSource(const Source &) {}

   virtual void write(std::string) {}


protected:

   XmlBuilder();

   xmlBase::XmlParser * m_parser;

   DOMDocument * m_doc;

};

} // namespace Likelihood

#endif // Likelihood_XmlBuilder_h
