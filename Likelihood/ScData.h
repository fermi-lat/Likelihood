/** 
 * @file ScData.h
 * @brief Declaration for ScData class, which contains the spacecraft data
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ScData.h,v 1.11 2003/10/22 16:31:33 jchiang Exp $
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
 * @brief Singleton container for ScNtuple data.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ScData.h,v 1.11 2003/10/22 16:31:33 jchiang Exp $
 */

class ScData {

public:

   ~ScData(){}

   /// Method to read in the spacecraft data.
   static void readData(std::string file, int hdu, bool clear = false);
   
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

   /// The spacecraft data itself. (This may be moved to the private 
   /// area and replaced here with access methods.)
   static std::vector<ScNtuple> vec;
#endif // SWIG

   /// Return the spacecraft z-axis as a function of MET.
   astro::SkyDir &zAxis(double time);

   /// Return the spacecraft x-axis as a function of MET.
   astro::SkyDir &xAxis(double time);

   /// Returns the Singleton object pointer.
   static ScData * instance();

protected:

   ScData(){}

private:

   static ScData * s_instance;
   
   static std::string s_scFile;
   static int s_scHdu;

   static double s_tstep;

   astro::SkyDir m_zAxis;
   astro::SkyDir m_xAxis;

};

} // namespace Likelihood
#endif // Likelihood_ScData_h
