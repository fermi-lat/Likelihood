/** 
 * @file LogLike.h
 * @brief Declaration of LogLike class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/LogLike.h,v 1.18 2005/03/02 04:51:10 jchiang Exp $
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/LogLike.h,v 1.18 2005/03/02 04:51:10 jchiang Exp $
 */

class LogLike : public SourceModel {
    
public:

   LogLike(const Observation & observation) : SourceModel(observation) {
      if (s_FT1_columns.size() == 0) {
         setFT1_columns();
      }
      deleteAllSources();
   }

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

   void getEvents(std::string event_file);

   void computeEventResponses(Source &src, double sr_radius = 30);

   void computeEventResponses(std::vector<DiffuseSource *> &srcs, 
                              double sr_radius = 30);

   void computeEventResponses(double sr_radius = 30);

   unsigned long nEvents() const {return m_events.size();}

   const std::vector<Event> & events() const {
      return m_events;
   }

protected:

   virtual LogLike * clone() const {
      return new LogLike(*this);
   }

   std::vector<Event> m_events;

private:

   Npred m_Npred;

   static std::vector<std::string> s_FT1_columns;

   static void setFT1_columns();

   void get_diffuse_names(tip::Table * events, 
                          std::vector<std::string> & names) const;

   double logSourceModel(const Event & event) const;

   void getLogSourceModelDerivs(const Event & event,
                                std::vector<double> & derivs) const;

};

} // namespace Likelihood

#endif // Likelihood_LogLike_h
