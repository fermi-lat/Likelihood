/**
 * @file BinnedLikelihood2.h
 * @brief Binned version of the log-likelihood function.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/BinnedLikelihood2.h,v 1.53 2011/04/01 21:16:47 jchiang Exp $
 */

#ifndef Likelihood_BinnedLikelihood2_h
#define Likelihood_BinnedLikelihood2_h

#include <map>

#include "tip/Image.h"

#include "optimizers/dArg.h"

#include "Likelihood/Accumulator.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/LogLike.h"
#include "Likelihood/Pixel.h"

namespace Likelihood {

   class Drm;
   class SourceMap;

/*
 * @class BinnedLikelihood2
 *
 */

class BinnedLikelihood2 : public LogLike {

public:

   BinnedLikelihood2(const CountsMap & dataMap, 
                    const Observation & observation,
                    const std::string & srcMapsFile="",
                    bool computePointSources=true,
                    bool applyPsfCorrections=true,
                    bool performConvolution=true,
                    bool resample=true,
                    double resamp_factor=2,
                    double minbinsz=0.1);

   BinnedLikelihood2(const BinnedLikelihood2 & other);

   virtual ~BinnedLikelihood2() throw();

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
                        bool addPointSources=true,
                        bool loadMaps=true,
                        bool createAllMaps=false);

   virtual CountsMap * createCountsMap() const;

   double npred();

   const SourceMap & sourceMap(const std::string & name) const {
      std::map<std::string, SourceMap *>::const_iterator srcMap
         = m_srcMaps.find(name);
      if (srcMap == m_srcMaps.end()) {
         throw std::runtime_error("Cannot find source map named: " + name);
      }
      return *(srcMap->second);
   }

   SourceMap * getSourceMap(const std::string & srcName,
                            bool verbose=true) const;

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

   /// The source component model-weighted counts spectrum, i.e.,
   /// based on the source model weights per pixel, this is the counts
   /// spectrum weighted by the probability of attributing counts in
   /// each pixel to this source.
   std::vector<double> countsSpectrum(const std::string & srcName) const;

   /// Predicited counts spectrum for the fixed model components
   /// summed together.
   std::vector<double> fixedModelSpectrum() const;

   virtual void addSource(Source * src);

   virtual Source * deleteSource(const std::string & srcName);

   virtual void syncParams();

   virtual void syncSrcParams(const std::string & srcName);

   virtual double NpredValue(const std::string & srcName) const;

   void set_klims(size_t kmin, size_t kmax) {
      m_modelIsCurrent = false;
      m_kmin = kmin;
      m_kmax = kmax;
      buildFixedModelWts();
   }

   std::pair<int, int> klims() const {
      return std::make_pair(static_cast<int>(m_kmin), static_cast<int>(m_kmax));
   }

   bool apply_drm() const {
      return m_apply_drm;
   }
   
   void set_apply_drm(bool apply_drm) {
      m_apply_drm = apply_drm;
      if (apply_drm && m_drm==0) {
         create_drm();
      }
   }
   
protected:

   virtual BinnedLikelihood2 * clone() const {
      return new BinnedLikelihood2(*this);
   }

   BinnedLikelihood2 & operator=(const BinnedLikelihood2 & rhs) {
      return *this;
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
   
   double m_minbinsz;

   std::vector<std::string> m_fixedSources;

   std::vector<std::pair<double, double> > m_fixedModelWts;
   std::map<std::string, double> m_fixedModelNpreds;

   /// Map of model parameters, to be used to determine if fixed
   /// sources have changed parameter values.
   std::map<std::string, std::vector<double> > m_modelPars;

   /// Accumulators for derivatives.
   mutable std::map<long, Accumulator> m_posDerivs;
   mutable std::map<long, Accumulator> m_negDerivs;

   /// Minimum and maximum energy plane indexes to use in likelihood 
   /// calculations.
   size_t m_kmin, m_kmax;

   /// Summed npred values at each energy boundary value for fixed sources.
   std::vector<double> m_fixedNpreds;

   /// Energy dispersion matrix
   Drm * m_drm;

   /// Flag to apply energy dispersion
   bool m_apply_drm;

   void createSourceMaps();

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

   void buildFixedModelWts();

   bool fixedModelUpdated() const;

   void create_drm();

};

}

#endif // Likelihood_BinnedLikelihood2_h
