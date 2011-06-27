/**
 * @file XmlParser.cxx
 * @brief Declaration of XmlParser::s_instance
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/XmlParser.cxx,v 1.3 2006/03/16 06:20:08 jchiang Exp $
 */

#define ST_DLL_EXPORTS
#include "Likelihood/XmlParser.h"
#undef ST_DLL_EXPORTS

namespace Likelihood {

xmlBase::XmlParser * XmlParser::s_instance(0);

xmlBase::XmlParser * XmlParser_instance() {
   return XmlParser::instance();
}

} // namespace Likelihood
