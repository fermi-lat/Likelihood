/**
 * @file EventContainer.h
 * @brief Container class for FT1 event data.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/EventContainer.h,v 1.4 2005/10/18 21:40:52 jchiang Exp $
 */

#ifndef Likelihood_EventContainer_h
#define Likelihood_EventContainer_h

#include <string>
#include <vector>

#include "Likelihood/Event.h"

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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/EventContainer.h,v 1.4 2005/10/18 21:40:52 jchiang Exp $
 */

class EventContainer {

public:

   EventContainer(const ResponseFunctions & respFuncs, 
                  const RoiCuts & roiCuts, const ScData & scData) 
      : m_respFuncs(respFuncs), m_roiCuts(roiCuts), m_scData(scData) {
      if (s_FT1_columns.size() == 0) {
         setFT1_columns();
      }
   }

   ~EventContainer() {}

   void getEvents(std::string event_file);
                  
   void computeEventResponses(Source & src, double sr_radius=80.);

   void computeEventResponses(std::vector<DiffuseSource *> &srcs, 
                              double sr_radius=80.);

   std::vector<double> nobs(const std::vector<double> & ebounds) const;

   unsigned long nEvents() const {
      return m_events.size();
   }

   const std::vector<Event> & events() const {
      return m_events;
   }

   void clear() {
      m_events.clear();
   }

private:
   
   const ResponseFunctions & m_respFuncs;
   const RoiCuts & m_roiCuts;
   const ScData & m_scData;

   std::vector<Event> m_events;

   static std::vector<std::string> s_FT1_columns;

   void setFT1_columns() const;

   void get_diffuse_names(tip::Table * events, 
                          std::vector<std::string> & names) const;

   std::string sourceName(const std::string & name) const;

};

} // namespace Likelihood

#endif // Likelihood_EventContainer_h
