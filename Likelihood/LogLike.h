/** 
 * @file LogLike.h
 * @brief Declaration of LogLike class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/LogLike.h,v 1.28 2006/09/18 20:59:24 jchiang Exp $
 */

#ifndef Likelihood_LogLike_h
#define Likelihood_LogLike_h

#include <map>

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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/LogLike.h,v 1.28 2006/09/18 20:59:24 jchiang Exp $
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

protected:

   virtual LogLike * clone() const {
      return new LogLike(*this);
   }

   mutable unsigned long m_nevals;

private:

   Npred m_Npred;

   std::map<std::string, double> m_npredValues;

   double logSourceModel(const Event & event) const;

   void getLogSourceModelDerivs(const Event & event,
                                std::vector<double> & derivs) const;

};

} // namespace Likelihood

#endif // Likelihood_LogLike_h
