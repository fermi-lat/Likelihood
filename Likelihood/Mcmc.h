/**
 * @file Mcmc.h
 * @brief Mcmc (Markov Chain Monte Carlo) class declaration
 * @author J. Chiang
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Mcmc.h,v 1.6 2003/06/11 17:08:02 jchiang Exp $
 */

#ifndef Likelihood_Mcmc_h
#define Likelihood_Mcmc_h

#include <vector>
#include <string>
#include "Likelihood/Parameter.h"
#include "Likelihood/Statistic.h"
#include "Likelihood/LikelihoodException.h"

namespace Likelihood {

/**
 * @class Mcmc
 *
 * @brief Apply the Variable-at-a-Time Metropolis-Hastings algorthim
 * to the (free) Parameters of a Statistic object.
 *
 * The transition probability distributions are specified along each
 * dimension by top-hat functions.  The widths of these top-hats may
 * either estimated, by default, as the rough 1-sigma error using an
 * approximate Hessian at the starting point via e.g., Minuit; or they
 * may be specified by hand through the setTransitionWidths method.
 * Because the Parameters are generally bounded, the transition
 * probabilities must be renormalized at each trial point by the
 * fraction of the top-hat contained within the boudaries, hence the
 * need for Metropolis-Hastings rather than simply the Metropolis
 * algorithm.
 *
 * Priors that are functions of the same set of Parameters in the form
 * of other Statistic objects can be applied.  As with the Statisic
 * itself, these priors need not be normalized to unity.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Mcmc.h,v 1.6 2003/06/11 17:08:02 jchiang Exp $
 */

class Mcmc {

public:

   Mcmc(Statistic &stat);
   ~Mcmc() {}

   void addPriors(std::vector<Statistic *> &priors) {m_priors = priors;}

   void generateSamples(std::vector< std::vector<double> > &samples,
                        unsigned long nsamp=1e4);

   //! Set the transition probablity widths by hand
   void setTransitionWidths(std::vector<double> &transitionWidths)
      {m_transitionWidths = transitionWidths;}

   //! Useful for restarting the MCMC
   void getTransitionWidths(std::vector<double> &transitionWidths)
      {transitionWidths = m_transitionWidths;}

   //! write samples to a FITS binary table
   static void writeSamples(std::string filename, 
                            std::vector< std::vector<double> > &samples)
      throw(LikelihoodException);

private:

   Statistic *m_stat;

   std::vector<Statistic *> m_priors;

   std::vector<double> m_transitionWidths;

   void estimateTransWidths();

   double drawValue(Parameter &param, double transitionWidth, 
                    double &transProbRatio);

};

} // namespace Likelihood

#endif // Likelihood_Mcmc_h
