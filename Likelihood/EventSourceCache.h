/** 
 * @file EventSourceCache.h
 * @brief Declaration of EventSourceCache class
 * @author S. Fegan
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/LogLike.h,v 1.37 2009/03/17 20:02:21 jchiang Exp $
 */

#ifndef Likelihood_EventSourceCache_h
#define Likelihood_EventSourceCache_h

#include<vector>
#include<map>

//#define DEBUG_EVENT_SOURCE_CACHE

namespace Likelihood {

template<typename T>
class EventSourceCache {

public:

  class EventRef
  {
  public:
    EventRef(): m_eventcache(), m_ievent() { }
    EventRef(EventSourceCache* eventcache, unsigned ievent): 
      m_eventcache(eventcache), m_ievent(ievent) { }
    inline T& getCachedValue(const std::string& srcname);
  private:
    EventSourceCache* m_eventcache;
    unsigned m_ievent;
  };

  EventSourceCache(unsigned nevents = 0): m_nevents(nevents), m_cache() { }

  ~EventSourceCache()
  {
#ifdef DEBUG_EVENT_SOURCE_CACHE
    std::cout << "Event/Source cache count:" << std::endl;
    for(std::map<std::string, unsigned>::iterator isrccnt=m_srcEvCount.begin();
	isrccnt!=m_srcEvCount.end(); isrccnt++)
      std::cout << isrccnt->first << ' ' << isrccnt->second << std::endl;
#endif
  }
    
  void clearAndResize(unsigned nevents = 0)
  {
    m_nevents = nevents;
    m_cache.clear();
  }

  void deleteSource(const std::string& srcname) { m_cache.erase(srcname); }

  EventRef getEventRef(unsigned ievent) { return EventRef(this,ievent); }

  T& getCachedValue(unsigned ievent,const std::string& srcname)
  {
    EvCache& evcache(m_cache[srcname]);
    if(evcache.size() != m_nevents)evcache.resize(m_nevents);
#ifdef DEBUG_EVENT_SOURCE_CACHE
    m_srcEvCount[srcname]++;
#endif
    return evcache.at(ievent);
  }

  std::vector<T>& getCachedEventValues(const std::string& srcname)
  {
    EvCache& evcache(m_cache[srcname]);
    if(evcache.size() != m_nevents)evcache.resize(m_nevents);
#ifdef DEBUG_EVENT_SOURCE_CACHE
    m_srcEvCount[srcname]+=m_nevents;
#endif
    return evcache;
  }

private:

  //EventSourceCache(const EventSourceCache&);
  EventSourceCache& operator= (const EventSourceCache&);

  typedef std::vector<T> EvCache;
  typedef std::map<std::string, EvCache> SrcEvCache;

  unsigned m_nevents;
  SrcEvCache m_cache;
#ifdef DEBUG_EVENT_SOURCE_CACHE
  std::map<std::string, unsigned> m_srcEvCount;
#endif
};

template<typename T> inline T& 
EventSourceCache<T>::EventRef::getCachedValue(const std::string& srcname)
{
  return m_eventcache->getCachedValue(m_ievent, srcname);
}

typedef std::pair<bool, double> CachedResponse;
typedef EventSourceCache<CachedResponse> ResponseCache;


} // namespace Likelihood

#endif // Likelihood_EventSourceCache_h
