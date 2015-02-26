/**
 * @file MultipleBrokenPowerLaw.h
 * @brief User configurable multiply broken power-law.
 * @author J. Chiang
 *
 * $Header: $
 */

#ifndef Likelihood_MultpleBrokenPowerLaw_h
#define Likelihood_MultpleBrokenPowerLaw_h

#include "optimizers/Arg.h"
#include "optimizers/Function.h"

namespace Likelihood {

/**
 * @class MultipleBrokenPowerLaw
 */

class MultipleBrokenPowerLaw : public optimizers::Function {

public:

   MultipleBrokenPowerLaw();

   double value(optimizers::Arg & x) const;
   
   double derivByParam(optimizers::Arg & x,
                       const std::string & paramName) const;

   double integral(optimizers::Arg & xmin, optimizers::Arg & xmax) const {
      // Temporary
      return 0;
   }

   virtual Function * clone() const {
      return new MultipleBrokenPowerLaw(*this);
   }

   void addParams(double norm, const std::vector<double> & photonIndexes,
                  const std::vector<double> & breakEnergies);

private:

   std::vector<std::string> m_indexNames;
   std::vector<double> m_breakEnergies;

   double norm(size_t k) const;

};

} // namespace Likelihood 

#endif // Likelihood_MultpleBrokenPowerLaw_h
