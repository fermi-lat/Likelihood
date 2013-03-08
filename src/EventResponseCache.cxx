#include "EventRespnseCache.h"

EventResponseCache::
EventResponseCache(Source::EDispMode edisp, unsigned nevent):
  m_resp(), m_resp_offset(), m_resp_count(), m_edisp(edisp) {
   switch(m_edisp) {
   case Source::ED_NONE:
      m_resp.resize(nevent);
      break;
   case Source::ED_GAUSSIAN:
      m_resp.resize(nevent*3);
      break;
   case Source::ED_FULL:
      m_resp_offset.resize(nevent);
      m_resp_count.resize(nevent);
      break;
   }
}

unsigned EventResponseCache::nevent() const {
   switch(m_edisp) {
   case Source::ED_NONE:
      return m_resp.size();
   case Source::ED_GAUSSIAN:
      return m_resp.size()/3;
   case Source::ED_FULL:
      return m_resp_offset.size();
   }
   assert(0);
}

void EventResponseCache::
getResponse(Source::Response& resp, unsigned ievent) const{
   switch(m_edisp) {
   case Source::ED_NONE:
      resp.resize(1);
      resp[0] = m_resp[ievent];
      break;
   case Source::ED_GAUSSIAN:
      resp.resize(3);
      std::copy(m_resp.begin() + ievent*3, m_resp.begin() + (ievent+1)*3,
		resp.begin());
      break;
   case Source::ED_FULL:
      resp.resize(m_resp_count[ievent]);
      std::copy(m_resp.begin() + m_resp_offset[ievent],
		m_resp.begin() + m_resp_offset[ievent] + m_resp_count[ievent],
		resp.begin());
      break;
   }
}

void EventResponseCache::
setResponse(const Source::Response& resp, unsigned ievent) {
   switch(m_edisp) {
   case Source::ED_NONE:
      ST_DEBUG_ASSERT(resp.size() == 1);
      m_resp[ievent] = resp[0];
      break;
   case Source::ED_GAUSSIAN:
      ST_DEBUG_ASSERT(resp.size() == 3);
      m_resp[ievent*3+0] = resp[0];
      m_resp[ievent*3+1] = resp[1];
      m_resp[ievent*3+2] = resp[2];
      break;
   case Source::ED_FULL:
      m_resp_offset[ievent] = m_resp.size();
      m_resp_count[ievent]  = resp.size();
      m_resp.insert(m_resp.end(), resp.begin(), resp.end());
      break;
   }
}

static EventResponseCache * EventResponseCache::
createFromCalculation(const Source * src, const std::vector<Event> & events,
		      unsigned nthreads) {
   unsigned nevent = events.size();
   EventResponseCache * resp_cache = 
     new EventResponseCache(src->edisp(), nevent);
   
   unsigned ievent;
#pragma omp parallel for private(ievent) num_threads(m_nthread) // schedule(dynamic,std::max(1,nevent/(50*m_nthread)))
   for(ievent=0; ievent<nevent; ievent++) {
      const Event & event = events[ievent];	
      EventResponse resp;
      src->computeResponse(resp, event);
#pragma omp critical
      resp_cache->setResponse(resp, ievent);
   }
   return resp_cache;
}

static EventResponseCache * EventResponseCache::
createFromDiffRespCols(const Source * src, 
		       const std::vector<Event> & events) {
   
}
