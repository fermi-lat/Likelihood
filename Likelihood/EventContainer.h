/**
 * @file EventContainer.h
 * @brief Container class for FT1 event data.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/EventContainer.h,v 1.9 2010/05/11 15:51:43 jchiang Exp $
 */

#ifndef Likelihood_EventContainer_h
#define Likelihood_EventContainer_h

#include <string>
#include <vector>

#include "Likelihood/DiffRespNames.h"
#include "Likelihood/Event.h"

namespace st_stream {
   class StreamFormatter;
}

namespace tip {
   class Table;
}

namespace Likelihood {
   
   class DiffuseSource;
   class ResponseFunctions;
   class RoiCuts;
   class ScData;
   class Source;

/**
 * @class EventContainer
 * @brief Container class for FT1 event data.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/EventContainer.h,v 1.9 2010/05/11 15:51:43 jchiang Exp $
 */

class EventContainer {

public:

   EventContainer(const ResponseFunctions & respFuncs, 
                  const RoiCuts & roiCuts, const ScData & scData);

   ~EventContainer();

   void getEvents(std::string event_file);
                  
   void computeEventResponses(Source & src, double sr_radius=30.);

   void computeEventResponses(std::vector<DiffuseSource *> &srcs, 
                              double sr_radius=30.);

   std::vector<double> nobs(const std::vector<double> & ebounds,
                            const Source * src=0) const;

   unsigned long nEvents() const {
      return m_events.size();
   }

   const std::vector<Event> & events() const {
      return m_events;
   }

   std::vector<Event> & events() {
      return m_events;
   }

   void clear() {
      m_events.clear();
   }

private:
   
   const ResponseFunctions & m_respFuncs;
   const RoiCuts & m_roiCuts;
   const ScData & m_scData;

   st_stream::StreamFormatter * m_formatter;

   std::vector<Event> m_events;

   static std::vector<std::string> s_FT1_columns;

   void setFT1_columns() const;

   void get_diffuse_names(tip::Table * events, 
                          std::vector<std::string> & names) const;

   void get_diffuse_names(tip::Table * events, 
                          DiffRespNames & diffRespNames) const;

   std::string sourceName(const std::string & name) const;

};

} // namespace Likelihood

#endif // Likelihood_EventContainer_h
