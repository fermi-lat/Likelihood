/**
 * @file XmlParser.cxx
 * @brief Declaration of XmlParser::s_instance
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/XmlParser.cxx,v 1.1 2005/11/16 20:00:33 jchiang Exp $
 */

#define ST_DLL_EXPORTS
#include "XmlParser.h"
#undef ST_DLL_EXPORTS

namespace Likelihood {

xmlBase::XmlParser * XmlParser::s_instance(0);

} // namespace Likelihood
