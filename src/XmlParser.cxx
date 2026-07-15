/**
 * @file XmlParser.cpp
 * @brief Declaration of XmlParser::s_instance
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/XmlParser.cxx,v 1.3 2006/03/16 06:20:08 jchiang Exp $
 */

#include "Likelihood/XmlParser.h"

namespace Likelihood {

xml_framework::SafeXmlParser* XmlParser::s_instance(nullptr);

xml_framework::SafeXmlParser* XmlParser_instance() {
   return XmlParser::instance();
}

} // namespace Likelihood
