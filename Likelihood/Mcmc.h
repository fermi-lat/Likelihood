/**
 * @file Mcmc.h
 * @brief Mcmc (Markov Chain Monte Carlo) class declaration
 * @author J. Chiang
 * $Header$
 */

#ifndef Mcmc_h
#define Mcmc_h

#include <vector>

namespace Likelihood {

class Statistic;

/**
 * @class Mcmc
 *
 * @brief Apply the Variable-at-a-Time Metropolis-Hastings algorthim
 * to the (free) Parameters of a Statistic object.
 *
 * The transition probability distributions are specified along each
 * dimension by 1D Gaussian functions.  The widths of these Gaussians
 * may either estimated, by default, using an approximate Hessian at
 * the starting point; or they may be specified by hand through the
 * setTransitionWidths method.  Because the Parameters are generally
 * bounded, the transition probabilities must be renormalized at each
 * trial point by the fraction of the Gaussian contained within the
 * boudaries, hence the need for Metropolis-Hastings rather than
 * simply the Metropolis algorithm.
 *
 * Priors that are functions of the same set of Parameters in the form
 * of other Statistic objects can be applied.  As with the Statisic
 * itself, these priors need not be normalized to unity.
 *
 * @author J. Chiang
 *
 * $Header$
 */

class Mcmc {

public:

   Mcmc(Statistic &stat) {m_stat = stat;}
   ~Mcmc() {}

   void addPriors(std::vector<Statistic> &priors) {m_priors = priors;}

   void setBurnIn(long burnIn) {m_burnIn = burnIn;}

   void generateSamples(std::vector< std::vector<double> > samples,
                        long nsamp=1e4);

   void setTransitionWidths(std::vector<double> transitionWidths)
      {m_transitionWidths = transitionWidths;}

private:

   Statistic m_stat;

   std::vector<Statistic> m_priors;

   long m_burnIn;

   std::vector<double> m_transitionWidths;
};

} // namespace Likelihood

#endif // Mcmc_h
