/** 
 * @file LogLike.h
 * @brief Declaration of LogLike class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/LogLike.h,v 1.36 2009/01/19 15:18:17 sfegan Exp $
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
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/LogLike.h,v 1.36 2009/01/19 15:18:17 sfegan Exp $
 */

class LogLike : public SourceModel {
    
public:

   LogLike(const Observation & observation);

   virtual ~LogLike() {}

   virtual double value(optimizers::Arg&) const;

   virtual double value() const {
      optimizers::Arg dummy;
      return value(dummy);
   }

   /// Return the derivatives wrt the free parameters, overloading
   /// the Function method
   virtual void getFreeDerivs(optimizers::Arg&, 
                              std::vector<double> &freeDerivs) const;

   virtual void getFreeDerivs(std::vector<double> &freeDerivs) const {
      optimizers::Arg dummy;
      getFreeDerivs(dummy, freeDerivs);
   }

   virtual void addSource(Source * src);

   virtual Source * deleteSource(const std::string & srcName);

   void getEvents(std::string event_file);

   void computeEventResponses(double sr_radius=30.);

   virtual void syncParams();

   virtual void syncSrcParams(const std::string & srcName);

   virtual double NpredValue(const std::string & srcName) const;

   void restoreBestFit();

   void saveCurrentFit();

protected:

   virtual LogLike * clone() const {
      return new LogLike(*this);
   }

   mutable unsigned long m_nevals;

   mutable double m_bestValueSoFar;

   void saveBestFit(double logLikeValue) const;

private:

   typedef Source::CachedResponse CachedResponse;

   Npred m_Npred;

   mutable Accumulator m_accumulator;

   std::map<std::string, double> m_npredValues;
   mutable std::vector<std::map<std::string, CachedResponse> > m_evSrcRespCache;

   double logSourceModel(const Event & event, 
			 std::map<std::string, CachedResponse>* srcRespCache=0)
     const;

   void getLogSourceModelDerivs(const Event & event,
                                std::vector<double> & derivs,
				std::map<std::string, CachedResponse>* srcRespCache=0) const;

   mutable std::vector<double> m_bestFitParsSoFar;

};

} // namespace Likelihood

#endif // Likelihood_LogLike_h
