/** 
 * @file LogLike.h
 * @brief Declaration of LogLike class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/LogLike.h,v 1.47 2016/03/29 23:53:49 echarles Exp $
 */

#ifndef Likelihood_LogLike_h
#define Likelihood_LogLike_h

#include <map>

#include "Likelihood/Accumulator.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/Event.h"
#include "Likelihood/Npred.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/SourceModel.h"
#include "Likelihood/EventSourceCache.h"

namespace tip {
   class Table;
}

namespace Likelihood {

/** 
 * @class LogLike
 *
 * @brief Objective function for the log(likelihood) of a model comprising
 * multiple Sources.
 *    
 */

class LogLike : public SourceModel {
    
public:

   LogLike(const Observation & observation);

   virtual ~LogLike() {}

   virtual double value(const optimizers::Arg&) const;

   virtual double value() const {
      optimizers::Arg dummy;
      return value(dummy);
   }

   /// Return the derivatives wrt the free parameters, overloading
   /// the Function method
   virtual void getFreeDerivs(const optimizers::Arg &,
                              std::vector<double> & freeDerivs) const;

   virtual void getFreeDerivs(std::vector<double> & freeDerivs) const {
      optimizers::Arg dummy;
      getFreeDerivs(dummy, freeDerivs);
   }

   virtual void addSource(Source * src, bool fromClone=true);

   virtual Source * deleteSource(const std::string & srcName);

   void getEvents(std::string event_file);

   void computeEventResponses(double sr_radius=30.);

   virtual void syncParams();

   virtual void syncSrcParams(const std::string & srcName);

   virtual double NpredValue(const std::string & srcName, bool weighted=false) const;

   void restoreBestFit();

   void saveCurrentFit();

   virtual void addPrior(size_t index, optimizers::Function & log_prior);

   virtual void set_ebounds(double emin, double emax);

   virtual void unset_ebounds();

protected:

   virtual LogLike * clone() const {
      return new LogLike(*this);
   }

   mutable unsigned long m_nevals;

   mutable double m_bestValueSoFar;

   bool m_use_ebounds;
   double m_emin;
   double m_emax;

   void saveBestFit(double logLikeValue) const;

protected:   

   mutable Kahan_Accumulator m_accumulator;

private:

   Npred m_Npred;

   std::map<std::string, double> m_npredValues;

   // Cache for instrument response to each event times source
   mutable ResponseCache m_respCache;

   double logSourceModel(const Event & event,
                         ResponseCache::EventRef* srcRespCache=0) const;

   void getLogSourceModelDerivs(const Event & event,
                                std::vector<double> & derivs,
                                ResponseCache::EventRef* srcRespCache=0) const;

   mutable std::vector<double> m_bestFitParsSoFar;

   void update_npreds();

};

} // namespace Likelihood

#endif // Likelihood_LogLike_h
