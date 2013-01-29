/**
 * @file BinnedLikelihood.h
 * @brief Binned version of the log-likelihood function.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/BinnedLikelihood.h,v 1.69 2013/01/28 12:40:38 sfegan Exp $
 */

#ifndef Likelihood_BinnedLikelihood_h
#define Likelihood_BinnedLikelihood_h

#include <map>
#include <stdexcept>

#include "tip/Image.h"

#include "optimizers/dArg.h"

#include "Likelihood/Accumulator.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/LogLike.h"
#include "Likelihood/Pixel.h"

namespace Likelihood {

class BinnedLikelihood : public LikelihoodBase {

public:

   BinnedLikelihood(const CountsMap & dataMap, 
                    const Observation & observation,
                    const std::string & srcMapsFile="",
                    bool computePointSources=true,
                    bool applyPsfCorrections=true,
                    bool performConvolution=true,
                    bool resample=true,
                    double resamp_factor=2,
                    double minbinsz=0.1);

   BinnedLikelihood(const BinnedLikelihood & other);

   virtual ~BinnedLikelihood() throw();

   /// Used by BinnedAnalysis
   const std::vector<double> & energies() const;
   const std::vector<double> & countsSpectrum() const;
   const std::vector<double> & modelCountsSpectrum() const;
   void set_klims(size_t kmin, size_t kmax);
   bool fixedModelUpdated() const;
   void buildFixedModelWts();
   const std::vector<double> & fixedModelSpectrum() const;
   const std::vector<double> & countsSpectrum(const std::string & srcName) const;

};

}

#endif // Likelihood_BinnedLikelihood_h
