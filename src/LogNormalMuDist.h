/**
 * @file LogNormalMuDist.h
 * @brief Function to compute a log-Normal distribution of sampling
 * points in cos(theta) for diffuse response integral.
 * 
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header$
 */

#ifndef Likelihood_LogNormalMuDist_h
#define Likelihood_LogNormalMuDist_h

#include <vector>

namespace Likelihood {

class LogNormalMuDist {

public:

   static LogNormalMuDist * instance();

   const std::vector<double> & muPoints(double energy) const;

protected:

   LogNormalMuDist(double muSlope=-0.4, double muIntercept=0.46,
                   double sigma=1.20, double emin=30, double emax=3e5,
                   size_t numEnergies=10);

private:

   double m_sigma;
   double m_logEmin;
   double m_logEmax;
   size_t m_numEnergies;
   double m_estep;

   std::vector< std::vector<double> > m_muPoints;

   static LogNormalMuDist * s_instance;

   void createSample(double mu, std::vector<double> & mus);

   size_t energyIndex(double energy) const;

};

} // namespace Likelihood

#endif // Likelihood_LogNormalMuDist_h
