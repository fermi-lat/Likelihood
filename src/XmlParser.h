/**
 * @file XmlParser.h
 * @brief Singleton wrapper for xmlBase::XmlParser instance
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/XmlParser.h,v 1.3 2006/03/16 00:35:57 jchiang Exp $
 */

#ifndef Likelihood_XmlParser_h
#define Likelihood_XmlParser_h

#include "xmlBase/XmlParser.h"

#include "st_facilities/libStApiExports.h"

namespace Likelihood {

/**
 * @class XmlParser
 * @brief Provides a Singleton wrapper a the static instance
 * of xmlBase::XmlParser.
 *
 * @author J. Chiang
 */

class SCIENCETOOLS_API XmlParser {

public:

   /// @return Pointer to the single xmlBase::XmlParser instance.
   static xmlBase::XmlParser * instance() {
      if (s_instance == 0) {
         s_instance = new xmlBase::XmlParser();
      }
      return s_instance;
   }

private:

   static xmlBase::XmlParser * s_instance;

};

// Opaque wrapper since linkage of exported symbols from windows dlls is
// all fouled up.
xmlBase::XmlParser * XmlParser_instance();

} // namespace Likelihood

#endif // Likelihood_XmlParser_h
