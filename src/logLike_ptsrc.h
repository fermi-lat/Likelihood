/** @file logLike_ptsrc.h
 * @brief Declaration of logLike_ptsrc class
 * $Header:
 */

#ifndef logLike_ptsrc_h
#define logLike_ptsrc_h

#include "../Likelihood/Event.h"
#include "../Likelihood/Statistic.h"
#include "../Likelihood/RoiCuts.h"
#include "../Likelihood/PointSource.h"

namespace Likelihood {

/** 
 * @class logLike_ptsrc
 *
 * @brief Objective function for the log(likelihood) of a 
 * single point source.
 *
 * @author J. Chiang
 *    
 * $Header: 
 */

class logLike_ptsrc : public Statistic {
    
public:

   logLike_ptsrc(){};
   virtual ~logLike_ptsrc(){};

   //! return the objective function value taking the free parameters 
   //! as the function argument
   double value(const std::vector<double> &paramVec);

   double evaluate_at(const Event &) const;

   void getEvents(const std::string &event_file, int hdu);

private:

   std::vector<Event> m_events;

};

} // namespace Likelihood

#endif // logLike_ptsrc_h
