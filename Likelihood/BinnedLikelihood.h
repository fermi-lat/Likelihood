/**
 * @file BinnedLikelihood.h
 * @brief Binned version of the log-likelihood function.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/BinnedLikelihood.h,v 1.13 2004/11/03 23:50:16 jchiang Exp $
 */

#ifndef Likelihood_BinnedLikelihood_h
#define Likelihood_BinnedLikelihood_h

#include <map>

#include "optimizers/dArg.h"

#include "Likelihood/CountsMap.h"
#include "Likelihood/LogLike.h"
#include "Likelihood/Pixel.h"

namespace tip {
   class Image;
}

namespace Likelihood {

class SourceMap;

/*
 * @class BinnedLikelihood
 *
 */

class BinnedLikelihood : public LogLike {

public:

   BinnedLikelihood(const CountsMap & dataMap, 
                    const std::string & srcMapsFile="");

//   BinnedLikelihood(const std::string & dataMapFile);
                 
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
      if (m_srcMapsFile == "") {
         createSourceMaps();
      } else {
         readSourceMaps();
      }
   }

   virtual CountsMap * createCountsMap() const;

   double npred();

   const SourceMap & sourceMap(const std::string & name) const {
      return *(m_srcMaps.find(name)->second);
   }

   void saveSourceMaps(const std::string & filename="");

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

   void readSourceMaps(std::string filename="");

   void computeModelMap(double & npred) const;

   void computeModelMap(std::vector<double> & modelMap) const;

   // Implement some rune-like tip arcana.
   void setImageDimensions(tip::Image * image, long * dims) const;

   void identifyFilledPixels();
   
   void fitsReportError(FILE *stream, int status) const;

   bool sourceMapExists(const std::string & srcName, 
                        const std::string & fitsFile) const;

   void replaceSourceMap(const std::string & srcName, 
                         const std::string & fitsFile) const;

   void addSourceMap(const std::string & srcName, 
                     const std::string & fitsFile) const;
};

}

#endif // Likelihood_BinnedLikelihood_h
