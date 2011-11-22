/**
 * @brief WcsMapLibrary.cxx
 * @brief Singleton class that contains the only copies of FITS
 * images.  It stores the FITS data in WcsMap2 objects and provides
 * access to these objects via pointers in an internal map keyed by
 * a string comprising the name of the FITS file and the extension.
 * @author J. Chiang 
 *
 * $Header$
 */

#include <utility>

#include "Likelihood/MapBase.h"
#include "Likelihood/WcsMap2.h"
#include "Likelihood/WcsMapLibrary.h"

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

WcsMap2 * WcsMapLibrary::wcsmap(const std::string & filename,
                                const std::string & extname) {
   std::string key(filename + "::" + extname);
   MapLibrary_t::const_iterator it(m_library.find(key));
   if (it != m_library.end()) {
      // Entry for map exists and it has not been deleted 
      // (via MapBase::deleteMap()).
      return it->second;
   }
   m_library[key] = new WcsMap2(filename, extname);
   return m_library[key];
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

void WcsMapLibrary::add_observer(const MapBase * observer) {
   m_observers.insert(std::make_pair(observer, 1));
}

void WcsMapLibrary::remove_observer(const MapBase * observer) {
   m_observers.erase(observer);
}

void WcsMapLibrary::notify() {
   std::map<const MapBase *, int>::iterator it(m_observers.begin());
   for ( ; it != m_observers.end(); ++it) {
      it->first->update();
   }
}

} // namespace Likelihood
