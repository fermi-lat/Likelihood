/**
 * @file MultipleBrokenPowerLaw.h
 * @brief User configurable multiply broken power-law.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/MultipleBrokenPowerLaw.h,v 1.2 2015/03/03 18:05:36 jchiang Exp $
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

   virtual Function * clone() const {
      return new MultipleBrokenPowerLaw(*this);
   }

   void addParams(double norm, const std::vector<double> & photonIndexes,
                  const std::vector<double> & breakEnergies);

protected:

   double value(const optimizers::Arg & x) const;
   
   double derivByParamImp(const optimizers::Arg & x,
                          const std::string & paramName) const;

private:

   std::vector<std::string> m_indexNames;
   std::vector<double> m_breakEnergies;

   double norm(size_t k) const;

};

} // namespace Likelihood 

#endif // Likelihood_MultpleBrokenPowerLaw_h
