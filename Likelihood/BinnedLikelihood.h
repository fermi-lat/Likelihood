/**
 * @file BinnedLikelihood.h
 * @brief Binned version of the log-likelihood function.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/BinnedLikelihood.h,v 1.75 2015/12/10 00:57:58 echarles Exp $
 */

#ifndef Likelihood_BinnedLikelihood_h
#define Likelihood_BinnedLikelihood_h

#include <map>
#include <stdexcept>

#include "tip/Image.h"

#include "optimizers/dArg.h"

#include "Likelihood/Accumulator.h"
#include "Likelihood/CountsMapBase.h"
#include "Likelihood/LogLike.h"
#include "Likelihood/Pixel.h"

namespace Likelihood {

   class Drm;
   class SourceMap;

/*
 * @class BinnedLikelihood
 * @brief Binned version of the log-Likelihood function.
 *
 */

class BinnedLikelihood : public LogLike {

public:

   BinnedLikelihood(CountsMapBase & dataMap, 
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

   virtual double value(optimizers::Arg &) const;

   virtual double value() const {
      optimizers::dArg dummy(0);
      return value(dummy);
   }

   virtual void getFreeDerivs(std::vector<double> & derivs) const;

   /// Create a counts map based on the current model.
   virtual CountsMapBase * createCountsMap(CountsMapBase & dataMap) const {
      std::vector<float> map;
      computeModelMap(map);
      dataMap.setImage(map);
      return &dataMap;
   }

   virtual void readXml(std::string xmlFile, 
                        optimizers::FunctionFactory & funcFactory,
                        bool requireExposure=true, 
                        bool addPointSources=true,
                        bool loadMaps=true);

   virtual CountsMapBase * createCountsMap() const;

   double npred();

   const SourceMap & sourceMap(const std::string & name) const;

   SourceMap & sourceMap(const std::string & name);

   SourceMap * getSourceMap(const std::string & srcName,
                            bool verbose=true) const;

   const CountsMapBase & countsMap() const {
      return m_dataMap;
   }

   void getNpreds(const std::string & srcName,
                  std::vector<double> & npreds) const;

   SourceMap * createSourceMap(const std::string & srcName);

   /// Instantiate or retrieve a SourceMap object for all sources and
   /// populate the internal source map hash with those objects.
   void createSourceMaps();

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
   std::vector<double> countsSpectrum(const std::string & srcName, 
                                      bool use_klims=true) const;

   /// Predicted counts spectrum for the fixed model components
   /// summed together.
   const std::vector<double> & fixedModelSpectrum() const;

   virtual void addSource(Source * src, bool fromClone=true);

   virtual Source * deleteSource(const std::string & srcName);
   
   void eraseSourceMap(const std::string & srcName);

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

   void setVerbose(bool verbose) {
      m_verbose = verbose;
   }

   void computeModelMap(std::vector<float> & modelMap) const;
   void updateModelMap(std::vector<float> & modeMap, 
                       const SourceMap * srcMap) const;
   void set_external_model_map(std::vector<float> * external_map) {
      m_external_model_map = external_map;
      if (external_map != 0) {
         external_map->resize(m_pixels.size()*(m_energies.size()-1));
      }
   }
   
   void setCountsMap(const std::vector<float> & counts);

   bool fixedModelUpdated() const;

   void buildFixedModelWts(bool process_all=false);

   void addFixedSource(const std::string & srcName);
   void deleteFixedSource(const std::string & srcName);

   const std::vector<double> & 
   modelCountsSpectrum(const std::string &srcname) const;

   void set_edisp_flag(bool use_edisp);

   bool use_edisp(const std::string & srcname="") const;

   void set_use_single_fixed_map(bool use_sfm);
   bool use_single_fixed_map() const;

   /// These are required since this inherits from LogLike rather than
   /// for SourceModel.  The inheritance hierarchy for this class and
   /// LogLike should be refactored.
   virtual void set_ebounds(double emin, double emax) {
      throw std::runtime_error("BinnedLikelihood::set_ebounds "
                               "not implemented.");
   }

   virtual void unset_ebounds() {
      throw std::runtime_error("BinnedLikelihood::unset_ebounds "
                               "not implemented.");
   }

protected:

   BinnedLikelihood & operator=(const BinnedLikelihood & rhs) {
      throw std::runtime_error("Copy-assignment operator not implemented");
   }

   virtual BinnedLikelihood * clone() const {
      return new BinnedLikelihood(*this);
   }

private:

   CountsMapBase& m_dataMap;

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

   bool m_verbose;

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

   bool m_use_edisp;
   Drm * m_drm;

   bool m_use_single_fixed_map;

   std::vector<float> * m_external_model_map;

   mutable std::map<std::string, std::vector<double> > m_true_counts;
   mutable std::map<std::string, std::vector<double> > m_meas_counts;

   mutable std::vector<double> m_fixed_counts_spec;

   std::map<std::string, std::map<size_t, size_t> > m_krefs;

   double computeModelMap() const;

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
   
   void replaceSourceMap_wcs(const SourceMap& srcMap, 
			     const std::string & fitsFile) const;

   void replaceSourceMap_healpix(const SourceMap& srcMap, 
				 const std::string & fitsFile) const;

   void appendSourceMap(const std::string & srcName, 
			const std::string & fitsFile) const;

   void appendSourceMap_wcs(const SourceMap& srcMap,
			    const std::string & fitsFile) const;

   void appendSourceMap_healpix(const SourceMap& srcMap, 
				const std::string & fitsFile) const;

   void computeCountsSpectrum();

   void computeCountsSpectrum_wcs();

   void computeCountsSpectrum_healpix();

   double spectrum(const Source * src, double energy) const;

   double pixelCounts(double emin, double emax, double y1, double y2) const;

   double NpredValue(const std::string & name, const SourceMap & srcMap) const;

   void computeFixedCountsSpectrum();

   void edisp_correction_factors(const std::string & srcName,
                                 const std::vector<double> & true_counts_spec,
                                 std::vector<double> &);

   Drm & drm();
};

}

#endif // Likelihood_BinnedLikelihood_h
