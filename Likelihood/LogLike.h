/** 
 * @file LogLike.h
 * @brief Declaration of LogLike class
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef Likelihood_LogLike_h
#define Likelihood_LogLike_h

#include "Likelihood/Event.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/logSrcModel.h"
#include "Likelihood/Npred.h"
#include "Likelihood/Table.h"

namespace Likelihood {

/** 
 * @class LogLike
 *
 * @brief Objective function for the log(likelihood) of a model comprising
 * multiple Sources.
 *
 * @author J. Chiang
 *    
 * $Header$
 */

class LogLike : public SourceModel {
    
public:

   LogLike() : m_eventData(0) {
      logSrcModel m_logSrcModel;
      Npred m_Npred;
      deleteAllSources();
   }

   virtual ~LogLike() {delete m_eventData;}

   virtual double value(optimizers::Arg&) const;

   /// Return the derivatives wrt the free parameters, overloading
   /// the Function method
   virtual void getFreeDerivs(optimizers::Arg&, 
                              std::vector<double> &freeDerivs) const;

   void getEvents(const std::string &event_file, int hdu);

//     void computeEventResponses(DiffuseSource &src, double sr_radius = 30);

//     void computeEventResponses(std::vector<DiffuseSource> &srcs, 
//                                double sr_radius = 30);

   void computeEventResponses(Source &src, double sr_radius = 30);

   void computeEventResponses(std::vector<DiffuseSource> &srcs, 
                              double sr_radius = 30);

   void computeEventResponses(double sr_radius = 30);

// Methods and data members from old Likelihood::Statistic:
   void readEventData(const std::string &eventFile, 
                      const std::string &colnames, int hdu);

   std::pair<long, double*> getEventColumn(const std::string &colname) const
      {return getColumn(*m_eventData, colname);}

   unsigned long nEvents() const {return m_events.size();}

protected:

   /// Generalized column access
   std::pair<long, double*> getColumn(const Table &tableData, 
                                      const std::string &colname) const
      throw(optimizers::Exception);

   /// Event data; read from m_eventFile, stored in Table form
   std::string m_eventFile;
   int m_eventHdu;
   Table *m_eventData;

private:

   std::vector<Event> m_events;

   // A bit of a kludge here making this guy mutable, but unavoidable
   // since getFreeDerivs is const in the base class and the
   // overloaded version needs to call the logSrcModel::mySyncParams
   // method.
   mutable logSrcModel m_logSrcModel;

   Npred m_Npred;

};

} // namespace Likelihood

#endif // Likelihood_LogLike_h
