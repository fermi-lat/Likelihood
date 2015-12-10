/**
 * @file WcsMapLibrary.h
 * @brief Singleton class that contains the only copies of FITS
 * images.  It stores the FITS data in ProjMap objects and provides
 * access to these objects via pointers in an internal map keyed by
 * the full path to the FITS file.
 * @author J. Chiang 
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/echarles/healpix_changes/Likelihood/Likelihood/WcsMapLibrary.h,v 1.4 2015/03/05 19:58:25 echarles Exp $
 */

#ifndef Likelihood_WcsMapLibrary_h
#define Likelihood_WcsMapLibrary_h

#include <map>
#include <string>

namespace Likelihood {

class MapBase;
// EAC, switch to using base class
class ProjMap;

class WcsMapLibrary {

public:

   ProjMap * wcsmap(const std::string & filename,
                    const std::string & extname);
   
   void delete_map(const std::string & filename,
                   const std::string & extname);

   bool has_map(const std::string & filename,
                const std::string & extname) const;

   void add_observer(MapBase * observer);

   void remove_observer(MapBase * observer);

   void notify();

   static WcsMapLibrary * instance() {
      if (s_instance == 0) {
         s_instance = new WcsMapLibrary();
      }
      return s_instance;
   }

   static void delete_instance() {
      delete s_instance;
      s_instance == 0;
   }

protected:

   WcsMapLibrary();

   ~WcsMapLibrary() throw();

private:

   typedef std::map<std::string, ProjMap *> MapLibrary_t;
   MapLibrary_t m_library;

   std::map<MapBase *, int> m_observers;

   static WcsMapLibrary * s_instance;
};

} // namespace Likelihood

#endif // Likelihood_WcsMapLibrary_h
