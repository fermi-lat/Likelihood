/** 
 * @file logLike_ptsrc.h
 * @brief Declaration of logLike_ptsrc class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/logLike_ptsrc.h,v 1.12 2003/08/06 20:52:09 jchiang Exp $
 */

#ifndef Likelihood_logLike_ptsrc_h
#define Likelihood_logLike_ptsrc_h

#include "Likelihood/Event.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/logSrcModel.h"
#include "Likelihood/Npred.h"
#include "Likelihood/Table.h"

namespace Likelihood {

/** 
 * @class logLike_ptsrc
 *
 * @brief Objective function for the log(likelihood) of a model comprising
 * multiple PointSources.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/logLike_ptsrc.h,v 1.12 2003/08/06 20:52:09 jchiang Exp $
 */

class logLike_ptsrc : public SourceModel {
    
public:

   logLike_ptsrc() : m_eventData(0) {
      logSrcModel m_logSrcModel;
      Npred m_Npred;
      deleteAllSources();
   }

   virtual ~logLike_ptsrc() {delete m_eventData;}

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

protected:

   //! generalized column access
   std::pair<long, double*> getColumn(const Table &tableData, 
                                      const std::string &colname) const
      throw(optimizers::Exception);

   //! Event data; read from m_eventFile, stored in Table form
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

#endif // Likelihood_logLike_ptsrc_h
