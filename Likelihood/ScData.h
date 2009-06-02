/** 
 * @file ScData.h
 * @brief Declaration for ScData class, which contains the spacecraft data
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ScData.h,v 1.26 2007/12/14 19:18:10 jchiang Exp $
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ScData.h,v 1.26 2007/12/14 19:18:10 jchiang Exp $
 */

class ScData {

public:

   ScData() {}

   ~ScData() {}

   /// Method to read in the spacecraft data.
   void readData(std::string file, double tstart, double tstop,
                 bool clear=false,
                 const std::string & sctable="SC_DATA");

   /// Read in data from several input files and check to see if no
   /// intervals are read in.
   void readData(const std::vector<std::string> & scFiles, 
                 double tstart, double tstop,
                 const std::string & sctable="SC_DATA");

   /// Return the spacecraft z-axis as a function of MET.
   astro::SkyDir zAxis(double time) const;

   /// Return the spacecraft x-axis as a function of MET.
   astro::SkyDir xAxis(double time) const;

   size_t numIntervals() const {
      return vec.size();
   }

   double start(size_t i) const {
      return vec.at(i).time;
   }

   double stop(size_t i) const {
      return vec.at(i).stoptime;
   }

   double livetime(size_t i) const {
      return vec.at(i).livetime;
   }

   const astro::SkyDir & zAxis(size_t i) const {
      return vec.at(i).zAxis;
   }

   unsigned int time_index(double time) const;

private:
   
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

   std::string m_scFile;

   static bool less_than_time(const ScNtuple & scDatum1,
                              const ScNtuple & scDatum2);

};

} // namespace Likelihood

#endif // Likelihood_ScData_h
