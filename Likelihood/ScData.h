/** 
 * @file ScData.h
 * @brief Declaration for ScData class, which contains the spacecraft data
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ScData.h,v 1.8 2003/07/19 04:38:02 jchiang Exp $
 */

#ifndef Likelihood_ScData_h
#define Likelihood_ScData_h

#include <string>
#include <vector>
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ScData.h,v 1.8 2003/07/19 04:38:02 jchiang Exp $
 */

class ScData {

public:

   ~ScData(){}

   //! method to read in the spacecraft data
   static void readData(const std::string &file, int hdu);
   
#ifndef SWIG
/** 
 * @class ScNtuple
 * @brief Nested NTuple class to represent spacecraft data.
 */
   class ScNtuple {
   public:
      ScNtuple(){}
      ~ScNtuple(){}
      double time;
      astro::SkyDir zenDir;
      astro::SkyDir xAxis;
      astro::SkyDir zAxis;
      int inSaa;
   };

   //! share the spacecraft data itself
   static std::vector<ScNtuple> vec;
#endif // SWIG

   //! returns the Singleton object pointer
   static ScData * instance();

protected:

   ScData(){}

private:

   static ScData * s_instance;
   
   static std::string s_scFile;
   static int s_scHdu;

};

} // namespace Likelihood
#endif // Likelihood_ScData_h
