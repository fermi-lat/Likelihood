/**
 * @file BinnedLikelihood2.h
 * @brief Binned version of the log-likelihood function.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/BinnedLikelihood2.h,v 1.9 2016/03/30 23:04:56 echarles Exp $
 */

#ifndef Likelihood_BinnedLikelihood2_h
#define Likelihood_BinnedLikelihood2_h

#include <map>

#include "tip/Image.h"

#include "optimizers/dArg.h"

#include "Likelihood/Accumulator.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/LogLike.h"

namespace Likelihood {

   class Drm;
   class SourceMap;

/*
 * @class BinnedLikelihood2
 *
 */

class BinnedLikelihood2 : public LogLike {

public:

   BinnedLikelihood2(const CountsMap & cmap,
                     const Observation & observation,
                     const std::string & smaps_file="",
                     bool computePointSources=true,
                     bool applyPsfCorrections=true,
                     bool performConvolution=true,
                     bool resample=true,
                     double resamp_factor=2,
                     double minbinsz=0.1);

   BinnedLikelihood2(const BinnedLikelihood2 & other);

   virtual ~BinnedLikelihood2() throw();

   virtual double value(const optimizers::Arg &) const;

   virtual double value() const {
      optimizers::dArg dummy(0);
      return value(dummy);
   }

   virtual void getFreeDerivs(std::vector<double> & derivs) const;

   virtual void readXml(std::string xmlFile, 
                        optimizers::FunctionFactory & funcFactory,
                        bool requireExposure=true, 
                        bool addPointSources=true,
                        bool loadMaps=true);

   virtual std::vector<double>::const_iterator setParamValues_(
      std::vector<double>::const_iterator);
   virtual std::vector<double>::const_iterator setFreeParamValues_(
      std::vector<double>::const_iterator);

   virtual void addSource(Source * src, bool fromClone=true);

   virtual Source * deleteSource(const std::string & srcName);

   virtual void syncParams();

   virtual void syncSrcParams(const std::string & srcName);

   virtual double NpredValue(const std::string & srcName, bool weighted = false) const;
   
   const std::vector<double> & energies() const {
      return m_energies;
   }

   const std::vector<double> & countsSpectrum() const {
      return m_countsSpectrum;
   }

   SourceMap * getSourceMap(const std::string & srcName,
                            bool verbose=true) const;

   SourceMap & sourceMap(const std::string & name) const {
      std::map<std::string, SourceMap *>::const_iterator srcMap
         = m_srcMaps.find(name);
      if (srcMap == m_srcMaps.end()) {
         throw std::runtime_error("Cannot find source map named: " + name);
      }
      return *(srcMap->second);
   }

   double npred() {
      return computeModelMap();
   }

   const CountsMap & countsMap() const {
      return m_cmap;
   }

   void getNpreds(const std::string & srcName,
                  std::vector<double> & npreds) const;

protected:

   virtual BinnedLikelihood2 * clone() const {
      return new BinnedLikelihood2(*this);
   }

   BinnedLikelihood2 & operator=(const BinnedLikelihood2 & rhs) {
      return *this;
   }

private:

   const CountsMap & m_cmap;
   std::string m_smaps_file;
   bool m_computePointSources;
   bool m_applyPsfCorrections;
   bool m_performConvolution;
   bool m_resample;
   double m_resamp_factor;
   double m_minbinsz;
   bool m_verbose;

   /// Energy dispersion matrix
   Drm * m_drm;

   mutable bool m_modelIsCurrent;

   std::vector<double> m_energies;
   std::vector<double> m_countsSpectrum;

   size_t m_npix, m_nee;

   mutable std::vector< std::vector<double> > m_model;

   std::map<std::string, SourceMap *> m_srcMaps;

   std::vector<std::string> m_fixedSources;

   typedef std::vector<std::vector<std::pair<double, double> > > ModelWeights_t;

   ModelWeights_t m_fixedModelWts;
   std::vector<double> m_fixedNpreds;
   std::map<std::string, double> m_fixedModelNpreds;
   std::map<std::string, std::vector<double> > m_modelPars;

   /// Accumulators for derivatives.
   mutable std::map<long, Accumulator> m_posDerivs;
   mutable std::map<long, Accumulator> m_negDerivs;

   void create_drm();

   double computeModelMap() const;

   void compute_model_wts(ModelWeights_t & modelWts);

   void addSourceWts(ModelWeights_t & modelWts,
                     const std::string & srcName,
		     SourceMap * srcMap=0,
                     bool subtract=false) const;

   double spectrum(const Source * src, double ee) const;

   double pixelCounts(double emin, double emax, double y1, double y2) const;

   bool fixedModelUpdated() const;

   void buildFixedModelWts();

   bool fileHasSourceMap(const std::string & srcName,
                         const std::string & fitsFile) const;

   double NpredValue(const std::string & name, SourceMap & srcMap) const;

   void computeCountsSpectrum();
};

} // namespace Likelihood

#endif // Likelihood_BinnedLikelihood2_h
