/**
 * @file CountsSpectra.h
 * @brief Encapsulation of counts spectra for a Likelihood fit.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/CountsSpectra.h,v 1.3 2006/09/14 20:38:41 peachey Exp $
 */

#ifndef Likelihood_CountsSpectra_h
#define Likelihood_CountsSpectra_h

#include <string>
#include <vector>

#include "Likelihood/BinnedLikelihood.h"
#include "Likelihood/BinnedLikelihood2.h"
#include "Likelihood/LogLike.h"

namespace Likelihood {

   class Observation;

/**
 * @class CountsSpectra
 *
 * @brief Class to handle the writing out of counts spectra for
 * unbinned and binned likelihood fits as a FITS binary table.
 *
 * @author J. Chiang
 */

class CountsSpectra {

public:

   CountsSpectra(const LogLike & logLike) 
      : m_logLike(logLike), m_observation(logLike.observation()), 
        m_binnedLike(0), m_binnedLike2(0) {
      m_binnedLike = 
         dynamic_cast<BinnedLikelihood *>(const_cast<LogLike *>(&m_logLike));
      if (m_binnedLike) {
         m_ebounds = m_binnedLike->energies();
         return;
      }
      m_binnedLike2 = 
         dynamic_cast<BinnedLikelihood2 *>(const_cast<LogLike *>(&m_logLike));
      if (m_binnedLike2) {
         m_ebounds = m_binnedLike2->energies();
      }
   }

   void setEbounds(double emin, double emax, size_t nbounds);

   void getSrcCounts(const std::string & srcName, 
                     std::vector<double> & srcCounts) const;

   void getSrcFluxes(const std::string & ptSrcName, 
                     std::vector<double> & srcFluxes) const;

   void getTotalSrcCounts(std::vector<double> & srcCounts) const;

   void getObsCounts(std::vector<double> & counts) const;

   void writeTable(const std::string & outfile) const;

   const std::vector<double> & ebounds() const {
      return m_ebounds;
   }

private:

   const LogLike & m_logLike;
   const Observation & m_observation;
   const BinnedLikelihood * m_binnedLike;
   const BinnedLikelihood2 * m_binnedLike2;
   std::vector<double> m_ebounds;

   void check_ebounds() const;

   void writeCounts(const std::string & outfile, 
                    const std::vector<std::string> & sourceNames) const;

   void writeFluxes(const std::string & outfile, 
                    const std::vector<std::string> & sourceNames) const;

   void writeEbounds(const std::string & outfile) const;

};

} // namespace Likelihood

#endif // Likelihood_CountsSpectra_h
