/** 
 * @file EventArg.h
 * @brief Declaration of EventArg class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/EventArg.h,v 1.6 2003/07/19 04:38:01 jchiang Exp $
 */

#ifndef Likelihood_EventArg_h
#define Likelihood_EventArg_h

#include "optimizers/Arg.h"
#include "Likelihood/Event.h"

namespace Likelihood {

/** 
 * @class EventArg
 *
 * @brief Concrete Arg subclass for encapsulating data of type Event.
 *
 * @authors J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/EventArg.h,v 1.6 2003/07/19 04:38:01 jchiang Exp $
 */

class EventArg : public optimizers::Arg {
    
public:
   
   EventArg(const Event &evt) : m_val(evt) {}
   virtual ~EventArg() {}

   void fetchValue(Event &evt) const {evt = m_val;}

private:

   Event m_val;

};

} // namespace Likelihood

#endif // Likelihood_EventArg_h
