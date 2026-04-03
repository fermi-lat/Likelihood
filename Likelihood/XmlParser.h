/**
 * @file XmlParser.h
 * @brief Singleton wrapper for SafeXmlParser instance
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/XmlParser.h,v 1.4 2006/03/16 06:20:08 jchiang Exp $
 */

#ifndef Likelihood_XmlParser_h
#define Likelihood_XmlParser_h

// RapidXML-based XML framework (replaces Xerces-C)
#include "xmlBase/safe_xml_parser.hpp"

#include "st_facilities/libStApiExports.h"

namespace Likelihood {

/**
 * @class XmlParser
 * @brief Provides a Singleton wrapper for the static instance
 * of SafeXmlParser.
 *
 * @author J. Chiang
 */
class SCIENCETOOLS_API XmlParser {

public:
   /// @return Pointer to the single SafeXmlParser instance.
   static xml_framework::SafeXmlParser* instance() {
      if (s_instance == nullptr) {
         s_instance = new xml_framework::SafeXmlParser();
      }
      return s_instance;
   }

   static void delete_instance() {
      delete s_instance;
      s_instance = nullptr;
   }

private:
   // Private constructor to prevent instantiation
   XmlParser() = default;
   ~XmlParser() = default;

   // Non-copyable
   XmlParser(const XmlParser&) = delete;
   XmlParser& operator=(const XmlParser&) = delete;

   static xml_framework::SafeXmlParser* s_instance;
};

// Opaque wrapper since linkage of exported symbols from windows dlls is
// all fouled up.
xml_framework::SafeXmlParser* XmlParser_instance();

} // namespace Likelihood

#endif // Likelihood_XmlParser_h
