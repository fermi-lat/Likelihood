/**
 * @file XmlBuilder.h
 * @brief Interface to builder classes for creating xml documents for
 * the flux and Likelihood packages.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/XmlBuilder.h,v 1.1 2004/02/20 00:02:03 jchiang Exp $
 */

#ifndef Likelihood_XmlBuilder_h
#define Likelihood_XmlBuilder_h

#include <xercesc/util/XercesDefs.hpp>

#include <string>

namespace xml {
   class XmlParser;
}

namespace Likelihood {

   class Source;

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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/XmlBuilder.h,v 1.1 2004/02/20 00:02:03 jchiang Exp $
 */

class XmlBuilder {

public:

   virtual ~XmlBuilder();

   virtual void addSource(Source &) {}

   virtual void write(std::string) {}

protected:

   XmlBuilder();

   xml::XmlParser * m_parser;

   XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument * m_doc;

};

} // namespace Likelihood

#endif // Likelihood_XmlBuilder_h
