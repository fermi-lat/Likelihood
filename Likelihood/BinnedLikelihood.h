/**
 * @file BinnedLikelihood.h
 * @brief Binned version of the log-likelihood function.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/BinnedLikelihood.h,v 1.2 2004/09/15 23:12:35 jchiang Exp $
 */

#ifndef Likelihood_BinnedLikelihood_h
#define Likelihood_BinnedLikelihood_h

#include "optimizers/dArg.h"

#include "Likelihood/CountsMap.h"
#include "Likelihood/Pixel.h"
#include "Likelihood/SourceModel.h"

namespace Likelihood {

/*
 * @class BinnedLikelihood
 *
 */

class BinnedLikelihood : public SourceModel {

public:

   BinnedLikelihood(const CountsMap & dataMap);
                 
   virtual ~BinnedLikelihood() throw() {}

   virtual double value(optimizers::Arg &) const;

   virtual double value() const {
      optimizers::dArg dummy(0);
      return value(dummy);
   }

   virtual void getFreeDerivs(std::vector<double> & derivs) const;

   virtual std::vector<double>::const_iterator setParamValues_(
      std::vector<double>::const_iterator);
   virtual std::vector<double>::const_iterator setFreeParamValues_(
      std::vector<double>::const_iterator);

private:

   const CountsMap & m_dataMap;

   std::vector<Pixel> m_pixels;
   std::vector<double> m_energies;

   mutable std::vector<double> m_model;
   mutable bool m_modelIsCurrent;

};

}

#endif // Likelihood_BinnedLikelihood_h
