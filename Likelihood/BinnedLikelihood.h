/**
 * @file BinnedLikelihood.h
 * @brief Binned version of the log-likelihood function.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/BinnedLikelihood.h,v 1.1 2004/09/13 15:30:39 jchiang Exp $
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
                 
   virtual ~BinnedLikelihood() throw() {
      try {
         delete m_modelMap;
      } catch (std::exception & eObj) {
         std::cerr << eObj.what() << std::endl;
      } catch (...) {
      }
   }

   virtual double value(optimizers::Arg &) const;

   virtual double value() const {
      optimizers::dArg dummy(0);
      return value(dummy);
   }

private:

   const CountsMap & m_dataMap;
   CountsMap * m_modelMap;

   std::vector<astro::SkyDir> m_pixelDirs;
   std::vector<double> m_pixelSolidAngles;
   std::vector<Pixel> m_pixels;
   std::vector<double> m_energies;

};

}

#endif // Likelihood_BinnedLikelihood_h
