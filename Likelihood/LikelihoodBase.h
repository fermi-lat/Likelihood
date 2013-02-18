/**
 * @file LikelihoodBase.h
 * @brief Base class for likelihood functions.
 * @author J. Chiang <jchiang@slac.stanford.edu>
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/Attic/LikelihoodBase.h,v 1.1.2.2 2013/01/30 16:09:29 jchiang Exp $
 */

#ifndef Likelihood_LikelihoodBase_h
#define Likelihood_LikelihoodBase_h

#include <map>
#include <string>
#include <vector>

#include <Accumulator.h>
#include "optimizers/Statistic.h"

namespace Likelihood {

class Observation;
class Source;

/**
 * @class LikelihoodBase
 */

class LikelihoodBase : public optimizers::Statistic {

   LikelihoodBase(const Observation & observation);

   virtual ~LikelihoodBase();
   
   double value() const;
   
   virtual void getFreeDerivs(std::vector<double> & derivs) const = 0;

   /// Member functions also needed by pyLikelihood classes.

   /// Used by AnalysisBase
   virtual void syncParams();

   virtual void setFreeParamValues(const std::vector<double> & pars);
   virtual void setParamValues(const std::vector<double> & pars);

   virtual void addSource(Source * src);
   virtual Source * getSource(const std::string & srcName);
   virtual Source * removeSource(const std::string & srcName);

   Source * deleteSource(const std::string & srcName) {
      // Would like to obsolete this function as its name suggests that the
      // source is deleted whereas it is simply removed from the analysis
      // and returned to the caller, who must delete it if required
      return removeSource(srcName);
   }

   /// @todo This member function should really be removed, but it is
   /// called by AnalysisBase.py.  AnalysisBase.py should simply call
   /// syncParams().
   virtual void syncSrcParams(const std::string & srcName) {
      (void)(srcName);
      syncParams();
   }

   virtual void getSrcNames(std::vector<std::string> & names) const;

   /// @todo This member function should be removed, but it is
   /// called by AnalysisBase.py.   Its purpose was to ensure that the
   /// log-likelihood for the initial set of parameters was saved.
   virtual void saveCurrentFit() {
      // do nothing.
   }

   virtual void restoreBestFit();

   virtual void writeXml(const std::string & xmlfile) const;
   
   virtual double NpredValue(const std::string & srcName) const = 0;

protected:

   virtual void addFreeSource(Source* src, const Source* callers_src);
   virtual void removeFreeSource(Source* src);

   /// Add source to the fixed source list and calculate contribution
   /// to all bins and to Npred
   virtual void addFixedSource(Source* src, const Source* callers_src);
   /// Remove source from fixed source list and trigger recalculation
   /// of fixed source contributions
   virtual void removeFixedSource(Source* src);
   /// Recalculate the cached fixed source values by clearing them all and
   /// repeatedly calling addFixedSources for all fixed sources
   virtual void rebuildFixedSourceCache();
   /// Test whether parametrs of fixed sources have changed (does not check
   /// the scale value)
   virtual bool hasFixedSourceModelChanged() const;

   /*
   class SourceInfo {
   public:
      Source *    source;
      size_t      parameter_ibegin;
      size_t      parameter_iend;
   };
   */

   const Observation &                      m_observation;

   std::map<std::string, Source *>          m_sources;

   std::vector<Source *>                    m_free_sources;
 
   std::map<Source *, std::vector<double> > m_fixed_source_param_values;
   std::vector<KahanAccumulator>            m_fixed_source_flux_in_bin;
   KahanAccumulator                         m_fixed_source_npred;
   KahanAccumulator                         m_zero_value;

   /// Best-fitting parameter values so far, as set by saveBestFit.
   mutable std::vector<double>              m_bestPars;

   /// Maximum log-likelihood value so far, as set by saveBestFit.
   mutable double                           m_maxValue;

   LikelihoodBase();

   LikelihoodBase(const LikelihoodBase & other) 
      : optimizers::Statistic(other) {}

   // This is called by LikelihoodBase::value().  Subclasses must
   // implement it to return the log-likelihood.
   double value(optimizers::Arg & ) const = 0;

   // These are pure virtual in optimizers::Function, so we
   // provide inaccessible concrete implementations.
   double derivByParam(optimizers::Arg &, const std::string &) const {
      return 0;
   }

   Function * clone() const {
      return 0;
   }

private:



};

} //namespace Likelihood

#endif // Likelihood_LikelihoodBase_h

