/** 
 * @file EventArg.h
 * @brief Declaration of EventArg class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/EventArg.h,v 1.5 2003/05/29 00:29:39 jchiang Exp $
 */

#ifndef Likelihood_EventArg_h
#define Likelihood_EventArg_h

#include "Likelihood/Arg.h"
#include "Likelihood/Event.h"

namespace Likelihood {

/** 
 * @class EventArg
 *
 * @brief Concrete Arg subclass for encapsulating data of type Event.
 *
 * @authors J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/EventArg.h,v 1.5 2003/05/29 00:29:39 jchiang Exp $
 */

class EventArg : public Arg {
    
public:
   
   EventArg(Event &evt) : m_val(evt) {}
   virtual ~EventArg() {}

   void fetchValue(Event &evt) const {evt = m_val;}

private:

   Event m_val;

};

} // namespace Likelihood

#endif // Likelihood_EventArg_h
