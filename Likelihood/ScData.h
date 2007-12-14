/** 
 * @file ScData.h
 * @brief Declaration for ScData class, which contains the spacecraft data
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ScData.h,v 1.25 2007/02/02 21:07:59 jchiang Exp $
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
 * @brief Container for spacecraft data.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ScData.h,v 1.25 2007/02/02 21:07:59 jchiang Exp $
 */

class ScData {

public:

   ScData() {}

   ~ScData() {}

   /// Method to read in the spacecraft data.
   void readData(std::string file, bool clear=false,
                 const std::string & sctable="SC_DATA");
   
   /// Method to read in the spacecraft data with start and stop times
   void readData(std::string file, double tstart, double tstop,
                 bool clear=false,
                 const std::string & sctable="SC_DATA");

   /// Read in data from several input files and check to see if no
   /// intervals are read in.
   void readData(const std::vector<std::string> & scFiles, 
                 double tstart, double tstop,
                 const std::string & sctable="SC_DATA");
   
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
      double stoptime;
      double livetime;
      astro::SkyDir zenDir;
      astro::SkyDir xAxis;
      astro::SkyDir zAxis;
      int inSaa;
   };

   /// The spacecraft data itself. (This may be moved to the private 
   /// area and replaced here with access methods.)
   std::vector<ScNtuple> vec;
#endif // SWIG

   /// Return the spacecraft z-axis as a function of MET.
   astro::SkyDir zAxis(double time) const;

   /// Return the spacecraft x-axis as a function of MET.
   astro::SkyDir xAxis(double time) const;

#ifndef SWIG
   /// Return a pair of iterators to the ScData intervals enclosing
   /// the desired start and end times.
   typedef std::vector<ScNtuple>::const_iterator Iterator;
   std::pair<Iterator, Iterator> bracketInterval(double startTime,
                                                 double stopTime) const;

   std::pair<Iterator, Iterator> 
   bracketInterval(const std::pair<double, double> & interval) const {
      return bracketInterval(interval.first, interval.second);
   }
#endif // SWIG

   unsigned int time_index(double time) const;

private:

   std::string m_scFile;
   int m_scHdu;

   static bool less_than_time(const ScNtuple & scDatum1,
                              const ScNtuple & scDatum2);

};

} // namespace Likelihood

#endif // Likelihood_ScData_h
