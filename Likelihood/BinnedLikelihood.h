/**
 * @file BinnedLikelihood.h
 * @brief Binned version of the log-likelihood function.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/BinnedLikelihood.h,v 1.6 2004/09/24 03:54:19 jchiang Exp $
 */

#ifndef Likelihood_BinnedLikelihood_h
#define Likelihood_BinnedLikelihood_h

#include <map>

#include "optimizers/dArg.h"

#include "Likelihood/CountsMap.h"
#include "Likelihood/Pixel.h"
#include "Likelihood/SourceModel.h"

namespace Likelihood {

class SourceMap;

/*
 * @class BinnedLikelihood
 *
 */

class BinnedLikelihood : public SourceModel {

public:

   BinnedLikelihood(const CountsMap & dataMap, 
                    const std::string & srcMapsFile="");
                 
   virtual ~BinnedLikelihood() throw() {}

   virtual double value(optimizers::Arg &) const;

   virtual double value() const {
      optimizers::dArg dummy(0);
      return value(dummy);
   }

   virtual void getFreeDerivs(std::vector<double> & derivs) const;

   /// Create a counts map based on the current model.
   virtual CountsMap * createCountsMap(const CountsMap & dataMap) const {
      return SourceModel::createCountsMap(dataMap);
   }

   virtual void readXml(std::string xmlFile, 
                        optimizers::FunctionFactory & funcFactory,
                        bool requireExposure=true) {
      SourceModel::readXml(xmlFile, funcFactory, requireExposure);
      createSourceMaps();
   }

   virtual CountsMap * createCountsMap() const;

   const SourceMap * sourceMap(const std::string & name) const {
      return m_srcMaps.find(name)->second;
   }

   void saveSourceMaps(std::string filename="") const;

   virtual std::vector<double>::const_iterator setParamValues_(
      std::vector<double>::const_iterator);
   virtual std::vector<double>::const_iterator setFreeParamValues_(
      std::vector<double>::const_iterator);

private:

   const CountsMap & m_dataMap;

   std::vector<Pixel> m_pixels;
   std::vector<double> m_energies;

   std::vector<unsigned int> m_filledPixels;

   std::map<std::string, SourceMap *> m_srcMaps;

   mutable std::vector<double> m_model;
   mutable bool m_modelIsCurrent;

   std::string m_srcMapsFile;

   void createSourceMaps();

   void computeModelMap(double & npred) const;

   void computeModelMap(std::vector<double> & modelMap) const;

   /// Implement some rune-like tip arcana.
   void setImageDimensions(tip::Image * image, long * dims) const;

   void identifyFilledPixels();

};

}

#endif // Likelihood_BinnedLikelihood_h
