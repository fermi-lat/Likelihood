/**
 * @file CountsSpectra.h
 * @brief Encapsulation of counts spectra for a Likelihood fit.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef Likelihood_CountsSpectra_h
#define Likelihood_CountsSpectra_h

#include <string>
#include <vector>

#include "Likelihood/BinnedLikelihood.h"
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
        m_binnedLike(0) {
      m_binnedLike = 
         dynamic_cast<BinnedLikelihood *>(const_cast<LogLike *>(&m_logLike));
      if (m_binnedLike) {
         m_ebounds = m_binnedLike->energies();
      }
   }

   void setEbounds(double emin, double emax, size_t nbounds);

   void getSrcCounts(const std::string & srcName, 
                     std::vector<double> & srcCounts) const;

   void getSrcFluxes(const std::string & ptSrcName, 
                     std::vector<double> & srcFluxes) const;

   void getObsCounts(std::vector<double> & counts) const;

   void writeTable(const std::string & outfile) const;

   const std::vector<double> & ebounds() const {
      return m_ebounds;
   }

private:

   const LogLike & m_logLike;
   const Observation & m_observation;
   const BinnedLikelihood * m_binnedLike;
   std::vector<double> m_ebounds;

   void check_ebounds() const;

};

} // namespace Likelihood

#endif // Likelihood_CountsSpectra_h
