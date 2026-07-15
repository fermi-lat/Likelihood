/**
 * @file XmlParser.h
 * @brief Singleton XML parser using RapidXML-based SafeXmlParser
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/XmlParser.h,v 1.3 2006/03/16 06:20:08 jchiang Exp $
 */

#ifndef Likelihood_XmlParser_h
#define Likelihood_XmlParser_h

#include "xmlBase/safe_xml_parser.hpp"

namespace Likelihood {

/**
 * @class XmlParser
 * @brief Singleton wrapper for SafeXmlParser
 */
class XmlParser {
public:
   static xml_framework::SafeXmlParser* instance() {
      if (s_instance == nullptr) {
         s_instance = new xml_framework::SafeXmlParser();
      }
      return s_instance;
   }

private:
   static xml_framework::SafeXmlParser* s_instance;
};

// Free function to get the parser instance (for C-style access)
xml_framework::SafeXmlParser* XmlParser_instance();

} // namespace Likelihood

#endif // Likelihood_XmlParser_h
