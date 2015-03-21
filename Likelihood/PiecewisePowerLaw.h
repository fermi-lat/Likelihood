/**
 * @file PiecewisePowerLaw.h
 * @brief User configurable piecewise continuous function made of
 * power-law segments.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/PiecewisePowerLaw.h,v 1.3 2015/03/16 03:26:57 jchiang Exp $
 */

#ifndef Likelihood_PiecewisePowerLaw_h
#define Likelihood_PiecewisePowerLaw_h

#include "optimizers/Arg.h"
#include "optimizers/Function.h"

namespace Likelihood {

/**
 * @class PiecewisePowerLaw
 */

class PiecewisePowerLaw : public optimizers::Function {

public:

   PiecewisePowerLaw();

   double integral(optimizers::Arg & xmin, optimizers::Arg & xmax) const {
      return 0;
   }

   virtual Function * clone() const {
      return new PiecewisePowerLaw(*this);
   }

   void addParams(double indexL, double indexH, 
                  const std::vector<double> & dNdEs,
                  const std::vector<double> & energies);

protected:

   double value(const optimizers::Arg & x) const;
   
   double derivByParamImp(const optimizers::Arg & x,
                          const std::string & paramName) const;

private:

   std::vector<std::string> m_dNdENames;
   std::vector<double> m_energies;
   
   bool m_decoupledNormPar;

   double plIndex(size_t k) const;
   double norm(size_t k) const;

   double derivByParam_decoupledNormPar(const optimizers::Arg & xarg,
                                        const std::string & paramName) const;

};

} // namespace Likelihood 

#endif // Likelihood_PiecewisePowerLaw_h
