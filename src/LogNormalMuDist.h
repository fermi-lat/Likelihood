/**
 * @file LogNormalMuDist.h
 * @brief Function to compute a log-Normal distribution of sampling
 * points in offset angle for diffuse response integral.
 * 
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/LogNormalMuDist.h,v 1.1 2007/12/11 07:14:49 jchiang Exp $
 */

#ifndef Likelihood_LogNormalMuDist_h
#define Likelihood_LogNormalMuDist_h

#include <vector>

namespace Likelihood {

/**
 * @class LogNormalDist
 */

class LogNormalMuDist {

public:

   static LogNormalMuDist * instance();

   const std::vector<double> & muPoints(double energy) const;

protected:

   LogNormalMuDist(double muSlope=-0.4, double muIntercept=0.46,
                   double sigma=1.20, double emin=30, double emax=3e5,
                   size_t numEnergies=10, size_t numMu=100);

private:

   double m_sigma;
   double m_logEmin;
   double m_logEmax;
   size_t m_numEnergies;
   size_t m_numMu;
   double m_estep;

   std::vector< std::vector<double> > m_muPoints;

   static LogNormalMuDist * s_instance;

   void createSample(double mu, std::vector<double> & mus);

   size_t energyIndex(double energy) const;

};

} // namespace Likelihood

#endif // Likelihood_LogNormalMuDist_h
