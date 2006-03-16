/**
 * @file XmlParser.cxx
 * @brief Declaration of XmlParser::s_instance
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/XmlParser.cxx,v 1.2 2006/03/13 17:14:16 jchiang Exp $
 */

#define ST_DLL_EXPORTS
#include "XmlParser.h"
#undef ST_DLL_EXPORTS

namespace Likelihood {

xmlBase::XmlParser * XmlParser::s_instance(0);

xmlBase::XmlParser * XmlParser_instance() {
   return XmlParser::instance();
}

} // namespace Likelihood
