/**
 * @file Mcmc.cxx
 * @brief Implementation for generating a Markov Chain Monte Carlo of
 * a Statistic object using the Variable-at-a-time Metropolis-Hastings
 * update method.
 * @author J. Chiang
 * $Header$
 */

#include <cmath>
#include "Likelihood/Mcmc.h"

namespace Likelihood {

Mcmc::Mcmc(Statistic &stat) {
   m_stat = &stat;
//   estimateTransWidths(stat, m_transitionWidths);
}

void Mcmc::generateSamples(std::vector< std::vector<double> > &samples, 
                           unsigned long nsamp) {
// Get initial values of Parameters... 
   std::vector<double> paramValues;
   m_stat->getFreeParamValues(paramValues);

// and the Parameter objects themselves (for the bounds information)
   std::vector<Parameter> params;
   m_stat->getFreeParams(params);

   samples.clear();

   while (samples.size() < nsamp) {
// Loop over parameters, treating each update step as a trial
      for (unsigned int i = 0; i < paramValues.size(); i++) {
         std::vector<double> newParamValues = paramValues;
         double transProbRatio;
         newParamValues[i] = drawValue(params[i], m_transitionWidths[i],
                                       transProbRatio);
// Hastings ratio
         double alpha = transProbRatio*exp(m_stat->value(newParamValues)
                                           - m_stat->value(paramValues));
// Metropolis rejection criterion (Use CLHEP call here?)
         double drand = static_cast<double>(rand())
            /static_cast<double>(RAND_MAX);
         if (drand < alpha) {
// Accept the new point in Parameter space
            paramValues[i] = newParamValues[i];
            params[i].setValue(paramValues[i]);
         }
// We always append the current point after the update step
         samples.push_back(paramValues);
      }
   }
}

void Mcmc::estimateTransWidths(Statistic *stat,
                               std::vector<double> &transitionWidths) {
// Use an approximate Hessian or call a Minuit routine to get estimates
// of each Parameters' error bars.
}

double Mcmc::drawValue(Parameter &param, double dx, double &transProbRatio) {
// Normalize the top-hat function (i.e., find the width) at the input
// Parameter location given the Parameter bounds.

   double xl = param.getBounds().first;
   double xu = param.getBounds().second;
   double x0 = param.getValue();
   double width = min(xu, x0 + dx) - max(xl, x0 - dx);

// Draw the trial value...
   double drand = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
   double y = drand*width + max(xl, x0 - dx);

// and compute the ratio of the transition probability densities
   transProbRatio = width/(min(xu, y + dx) - max(xl, y - dx));

   return y;
}

} // namespace Likelihood
