/** @file logLike_ptsrc.h
 * @brief Declaration of logLike_ptsrc class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/logLike_ptsrc.h,v 1.8 2003/04/25 18:32:20 jchiang Exp $
 */

#ifndef logLike_ptsrc_h
#define logLike_ptsrc_h

#include "Likelihood/Event.h"
#include "Likelihood/Statistic.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/logSrcModel.h"
#include "Likelihood/Npred.h"

namespace Likelihood {

/** 
 * @class logLike_ptsrc
 *
 * @brief Objective function for the log(likelihood) of a model comprising
 * multiple PointSources.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/logLike_ptsrc.h,v 1.8 2003/04/25 18:32:20 jchiang Exp $
 */

class logLike_ptsrc : public Statistic {
    
public:

   logLike_ptsrc() {
      logSrcModel m_logSrcModel;
      Npred m_Npred;
      deleteAllSources();
   }
   virtual ~logLike_ptsrc() {}

   //! return the objective function value taking the free parameters 
   //! as the function argument
   double value(const std::vector<double> &paramVec);

   //! return the derivatives wrt the free parameters
   void getFreeDerivs(std::vector<double> &freeDerivs);

   void getEvents(const std::string &event_file, int hdu);

//     void computeEventResponses(DiffuseSource &src, double sr_radius = 30);

//     void computeEventResponses(std::vector<DiffuseSource> &srcs, 
//                                double sr_radius = 30);

   void computeEventResponses(Source &src, double sr_radius = 30);

   void computeEventResponses(std::vector<DiffuseSource> &srcs, 
                              double sr_radius = 30);

   void computeEventResponses(double sr_radius = 30);

private:

   std::vector<Event> m_events;

   logSrcModel m_logSrcModel;

   Npred m_Npred;

};

} // namespace Likelihood

#endif // logLike_ptsrc_h
