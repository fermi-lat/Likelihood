/** 
 * @file EventArg.h
 * @brief Declaration of EventArg class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/EventArg.h,v 1.4 2003/03/17 00:53:42 jchiang Exp $
 */

#ifndef EventArg_h
#define EventArg_h

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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/EventArg.h,v 1.4 2003/03/17 00:53:42 jchiang Exp $
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

#endif // EventArg_h
