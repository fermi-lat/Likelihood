#ifndef ScData_h
#define ScData_h

//#include "../Likelihood/ScNtuple.h"
#include "astro/SkyDir.h"

namespace Likelihood {

/** 
 * @class ScData
 *
 * @brief Container for ScNtuple data.  Singleton.  Used by Response
 * and Source hierarchies.
 *
 * @author J. Chiang
 *    
 * $Header: */

class ScData {

public:

   ~ScData(){};

   //! method to read in the spacecraft data
   static void readData(const std::string &file, int hdu);
   
/** 
 * @class ScNtuple
 * @brief Nested NTuple class to represent spacecraft data.
 */
   class ScNtuple {
   public:
      ScNtuple(){};
      ~ScNtuple(){};
      double time;
      astro::SkyDir zenDir;
      astro::SkyDir xAxis;
      astro::SkyDir zAxis;
      int inSaa;
   };

   //! share the spacecraft data itself
   static std::vector<ScNtuple> vec;

   //! returns the Singleton object pointer
   static ScData * instance();

protected:

   ScData(){};

private:

   static ScData * s_instance;
   
   static std::string m_scFile;
   static int m_scHdu;

};

} // namespace Likelihood
#endif // ScData_h
