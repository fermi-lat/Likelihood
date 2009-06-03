/** 
 * @file ScData.h
 * @brief Spacecraft data class
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ScData.h,v 1.28 2009/06/02 20:51:39 jchiang Exp $
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
 * @brief Interface to spacecraft data
 *
 */

class ScData {

public:

   ScData() {}

   ~ScData() {}

   /// Method to read in the spacecraft data.
   void readData(std::string file, double tstart, double tstop,
                 bool clear=false, const std::string & sctable="SC_DATA");

   /// Read in data from several input files and check to see if no
   /// intervals are read in.
   void readData(const std::vector<std::string> & scFiles, 
                 double tstart, double tstop,
                 const std::string & sctable="SC_DATA");

   /// Livetime fraction as a function of MET
   double livetimefrac(double time) const;

   /// Spacecraft z-axis as a function of MET.
   astro::SkyDir zAxis(double time) const;

   /// Spacecraft x-axis as a function of MET.
   astro::SkyDir xAxis(double time) const;

   size_t numIntervals() const {
      return m_start.size();
   }

   double start(size_t i) const {
      return m_start.at(i);
   }

   double stop(size_t i) const {
      return m_stop.at(i);
   }

   double livetime(size_t i) const {
      return m_livetime.at(i);
   }

   const astro::SkyDir & xAxis(size_t i) const {
      return m_xAxis.at(i);
   }

   const astro::SkyDir & zAxis(size_t i) const {
      return m_zAxis.at(i);
   }

   size_t time_index(double time) const;

private:

   std::vector<double> m_start;
   std::vector<double> m_stop;
   std::vector<double> m_livetime;
   std::vector<astro::SkyDir> m_zAxis;
   std::vector<astro::SkyDir> m_xAxis;

   void clear_arrays();

};

} // namespace Likelihood

#endif // Likelihood_ScData_h
