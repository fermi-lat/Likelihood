/**
 * @brief WcsMapLibrary.cxx
 * @brief Singleton class that contains the only copies of FITS
 * images.  It stores the FITS data in ProjMap objects and provides
 * access to these objects via pointers in an internal map keyed by
 * a string comprising the name of the FITS file and the extension.
 * @author J. Chiang 
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/echarles/healpix_changes/Likelihood/src/WcsMapLibrary.cxx,v 1.3 2015/03/03 06:00:01 echarles Exp $
 */

#include <utility>

#include "Likelihood/MapBase.h"
#include "Likelihood/WcsMap2.h"
#include "Likelihood/HealpixProjMap.h"
#include "Likelihood/WcsMapLibrary.h"
#include "Likelihood/AppHelpers.h"

namespace Likelihood {

WcsMapLibrary * WcsMapLibrary::s_instance(0);

WcsMapLibrary::WcsMapLibrary() {}

WcsMapLibrary::~WcsMapLibrary() throw() {
   try {
      MapLibrary_t::iterator it(m_library.begin());
      for ( ; it != m_library.end(); ++it) {
         delete it->second;
      }
   } catch(...) {
   }
}

ProjMap * WcsMapLibrary::wcsmap(const std::string & filename,
                                const std::string & extname) {
   std::string key(filename + "::" + extname);
   MapLibrary_t::const_iterator it(m_library.find(key));
   if (it != m_library.end()) {
      return it->second;
   }
   astro::ProjBase::Method method = AppHelpers::checkProjectionMethod(filename,extname);
   ProjMap* theMap(0);
   switch ( method ) {
   case astro::ProjBase::WCS:
     theMap = new WcsMap2(filename, extname);
     break;
   case astro::ProjBase::HEALPIX:
     theMap = new HealpixProjMap(filename, extname);     
     break;
   default:
     break;
   }
   if ( theMap == 0 ) {
     std::string errMsg("Did not recognize ProjMap type at: ");
     errMsg += filename;
     if ( extname.size() > 0 ) {
       errMsg += "[" + extname + "]";
     }
     throw std::runtime_error(errMsg);
   }
   m_library[key] = theMap;
   return theMap;
}

void WcsMapLibrary::delete_map(const std::string & filename,
                               const std::string & extname) {
   std::string key(filename + "::" + extname);
   MapLibrary_t::const_iterator it(m_library.find(key));
   if (it != m_library.end()) {
      delete it->second;
      m_library.erase(key);
   }
   notify();
}

bool WcsMapLibrary::has_map(const std::string & filename,
                            const std::string & extname) const {
   std::string key(filename + "::" + extname);
   return m_library.find(key) != m_library.end();
}

void WcsMapLibrary::add_observer(MapBase * observer) {
   m_observers.insert(std::make_pair(observer, 1));
}

void WcsMapLibrary::remove_observer(MapBase * observer) {
   m_observers.erase(observer);
}

void WcsMapLibrary::notify() {
   std::map<MapBase *, int>::iterator it(m_observers.begin());
   for ( ; it != m_observers.end(); ++it) {
      it->first->update();
   }
}

} // namespace Likelihood
