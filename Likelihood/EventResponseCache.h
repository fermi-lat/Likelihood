/**
 * @file EventResponseCache.h
 * @brief Store source event responses for UnbinnedLikelihood
 * @author S. Fegan <sfegan@llr.in2p3.fr>
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/Attic/UnbinnedLikelihood.h,v 1.1.2.3 2013/03/06 14:36:50 sfegan Exp $
 */

#ifndef Likelihood_EventResponseCache_h
#define Likelihood_EventResponseCache_h

#include <vector>
#include <string>

#include "Event.h"
#include "Source.h"

namespace Likelihood {

class EventResponseCache {
public:
   typedef Source::Response Response;

   EventResponseCache(Source::EDispMode edisp, unsigned nevent);
   unsigned nevent() const;
   unsigned nbytes() const { return m_resp.size()*sizeof(double); }
   void getResponse(Response& resp, unsigned ievent) const;
   void setResponse(const Response& resp, unsigned ievent);
   
   static EventResponseCache * 
     createFromCalculation(const Source * src, 
			   const std::vector<Event> & events,
			   unsigned nthreads = 1);

   static EventResponseCache * 
     createFromDiffRespCols(const Source * src, 
			    const std::vector<Event> & events);

protected:
   std::vector<double>   m_resp;
   std::vector<unsigned> m_resp_offset;
   std::vector<unsigned> m_resp_count;
   Source::EDispMode     m_edisp;
};

} // namespace Likelihood

#endif // Likelihood_EventResponseCache_h
