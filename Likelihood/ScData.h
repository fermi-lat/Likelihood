/** 
 * @file ScData.h
 * @brief Declaration for ScData class, which contains the spacecraft data
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ScData.h,v 1.15 2005/01/23 00:38:07 jchiang Exp $
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ScData.h,v 1.15 2005/01/23 00:38:07 jchiang Exp $
 */

class ScData {

public:

   ~ScData(){}

   /// Method to read in the spacecraft data.
   static void readData(std::string file, bool clear=false);
   
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
   const astro::SkyDir &zAxis(double time);

   /// Return the spacecraft x-axis as a function of MET.
   const astro::SkyDir &xAxis(double time);

   /// Returns the Singleton object pointer.
   static ScData * instance();

#ifndef SWIG
   /// Return a pair of iterators to the ScData intervals enclosing
   /// the desired start and end times.
   typedef std::vector<ScNtuple>::iterator Iterator;
   static std::pair<Iterator, Iterator> bracketInterval(double startTime,
                                                        double stopTime);
   static std::pair<Iterator, Iterator> 
   bracketInterval(const std::pair<double, double> & interval) {
      return bracketInterval(interval.first, interval.second);
   }
#endif // SWIG

protected:

   ScData(){}

private:

   static ScData * s_instance;
   
   static std::string s_scFile;
   static int s_scHdu;

   static double s_tstep;

   astro::SkyDir m_zAxis;
   astro::SkyDir m_xAxis;

   static bool less_than_time(const ScNtuple & scDatum1,
                              const ScNtuple & scDatum2);

};

} // namespace Likelihood
#endif // Likelihood_ScData_h
