/**
 * @file BinnedLikelihood2.h
 * @brief Binned version of the log-likelihood function.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/BinnedLikelihood2.h,v 1.1 2011/06/27 05:05:37 jchiang Exp $
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

   BinnedLikelihood2(const CountsMap & cmap;
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

   virtual double value(optimizers::Arg &) const;

   virtual void getFreeDerivs(std::vector<double> & derivs) const;

   virtual void readXml(std::string xmlFile, 
                        optimizers::FunctionFactory & funcFactory,
                        bool requireExposure=true, 
                        bool addPointSources=true,
                        bool loadMaps=true,
                        bool createAllMaps=false);

   virtual std::vector<double>::const_iterator setParamValues_(
      std::vector<double>::const_iterator);
   virtual std::vector<double>::const_iterator setFreeParamValues_(
      std::vector<double>::const_iterator);

   virtual void addSource(Source * src);

   virtual Source * deleteSource(const std::string & srcName);

   virtual void syncParams();

   virtual void syncSrcParams(const std::string & srcName);

   virtual double NpredValue(const std::string & srcName) const;
   
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

   /// Energy dispersion matrix
   Drm * m_drm;

   bool m_modelIsCurrent;

   std::vector<double> m_energies;

   size_t m_npix, m_nee;

   std::vector< std::vector<double> > m_model;

   std::map<std::string, SourceMap *> m_srcMaps;

   std::vector<std::string> m_fixedSources;

   typedef std::vector<std::vector<std::pair<double, double> > > ModelWeights_t;

   ModelWeights_t m_fixedModelWts;

   /// Accumulators for derivatives.
   mutable std::map<long, Accumulator> m_posDerivs;
   mutable std::map<long, Accumulator> m_negDerivs;

   void create_drm();

   double computeModelMap() const;

   void update_model_wts();

   void addSourceWts(ModelWeights_t & modelWts,
                     const std::string & srcName,
                     const SourceMap * srcMap=0,
                     bool subtract=false) const;

   double spectrum(Sourc * src, double ee) const;

   double pixelCounts(double emin, double emax, double y1, double y2) const;
};

}

#endif // Likelihood_BinnedLikelihood2_h
