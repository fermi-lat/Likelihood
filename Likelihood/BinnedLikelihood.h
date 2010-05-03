/**
 * @file BinnedLikelihood.h
 * @brief Binned version of the log-likelihood function.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/BinnedLikelihood.h,v 1.39 2010/05/03 18:25:12 jchiang Exp $
 */

#ifndef Likelihood_BinnedLikelihood_h
#define Likelihood_BinnedLikelihood_h

#include <map>

#include "tip/Image.h"

#include "optimizers/dArg.h"

#include "Likelihood/CountsMap.h"
#include "Likelihood/LogLike.h"
#include "Likelihood/Pixel.h"

namespace Likelihood {

   class SourceMap;

/*
 * @class BinnedLikelihood
 * @brief Binned version of the log-Likelihood function.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/BinnedLikelihood.h,v 1.39 2010/05/03 18:25:12 jchiang Exp $
 */

class BinnedLikelihood : public LogLike {

public:

   BinnedLikelihood(const CountsMap & dataMap, 
                    const Observation & observation,
                    const std::string & srcMapsFile="",
                    bool computePointSources=true,
                    bool applyPsfCorrections=true,
                    bool performConvolution=true,
                    bool resample=true,
                    double resamp_factor=2);

   virtual ~BinnedLikelihood() throw();

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
                        bool requireExposure=true, 
                        bool addPointSources=true);

//   virtual CountsMap * createCountsMap() const;

   double npred();

   const SourceMap & sourceMap(const std::string & name) const {
      std::map<std::string, SourceMap *>::const_iterator srcMap
         = m_srcMaps.find(name);
      if (srcMap == m_srcMaps.end()) {
         throw std::runtime_error("Cannot find source map named: " + name);
      }
      return *(srcMap->second);
   }

   const CountsMap & countsMap() const {
      return m_dataMap;
   }

   void getNpreds(const std::string & srcName,
                  std::vector<double> & npreds) const;

   SourceMap * createSourceMap(const std::string & srcName);

   void saveSourceMaps(const std::string & filename="");

   virtual std::vector<double>::const_iterator setParamValues_(
      std::vector<double>::const_iterator);
   virtual std::vector<double>::const_iterator setFreeParamValues_(
      std::vector<double>::const_iterator);

   const std::vector<double> & energies() const {
      return m_energies;
   }

   const std::vector<double> & countsSpectrum() const {
      return m_countsSpectrum;
   }

   virtual void addSource(Source * src);

   virtual Source * deleteSource(const std::string & srcName);

   virtual void syncParams();

   virtual void syncSrcParams(const std::string & srcName);

   virtual double NpredValue(const std::string & srcName) const;

protected:

   virtual BinnedLikelihood * clone() const {
      return new BinnedLikelihood(*this);
   }

private:

   const CountsMap & m_dataMap;

   const std::vector<Pixel> & m_pixels;
   std::vector<double> m_energies;

   std::vector<double> m_countsSpectrum;

   std::vector<unsigned int> m_filledPixels;

   std::map<std::string, SourceMap *> m_srcMaps;

   mutable std::vector<double> m_model;
   mutable bool m_modelIsCurrent;

   std::string m_srcMapsFile;

   bool m_computePointSources;

   bool m_applyPsfCorrections;

   bool m_performConvolution;

   bool m_resample;
   
   double m_resamp_factor;

   std::vector<std::string> m_fixedSources;

   std::vector<std::pair<double, double> > m_fixedModelWts;
   std::map<std::string, double> m_fixedModelNpreds;
   
   void createSourceMaps();

   void readSourceMaps(std::string filename="");

   SourceMap * getSourceMap(const std::string & srcName) const;

   void computeModelMap(double & npred) const;

   void computeModelMap(std::vector<float> & modelMap) const;

   void addSourceWts(std::vector<std::pair<double, double> > & modelWts,
                     const std::string & srcName,
                     const SourceMap * srcMap=0, 
                     bool subtract=false) const;

   void setImageDimensions(tip::Image * image, long * dims) const;

   void identifyFilledPixels();
   
   bool fileHasSourceMap(const std::string & srcName, 
                         const std::string & fitsFile) const;

   void replaceSourceMap(const std::string & srcName, 
                         const std::string & fitsFile) const;

   void appendSourceMap(const std::string & srcName, 
                        const std::string & fitsFile) const;

   void computeCountsSpectrum();

   double spectrum(const Source * src, double energy) const;

   double pixelCounts(double emin, double emax, double y1, double y2) const;

   double NpredValue(const std::string & name, const SourceMap & srcMap) const;

   void updateFixedModelWts();
};

}

#endif // Likelihood_BinnedLikelihood_h
