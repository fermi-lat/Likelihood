/** 
 * @file ResponseFunctions.h
 * @brief A singleton class to contain the instrument response functions.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Response.h,v 1.21 2003/07/21 22:14:56 jchiang Exp $
 */

#ifndef Likelihood_ResponseFunctions_h
#define Likelihood_ResponseFunctions_h

#include <map>

namespace Likelihood {

/** 
 * @class ResponseFunctions
 *
 * @brief This class provides global access to a map of pointers to
 * latResponse::Irfs objects.  These pointers are indexed by event
 * type, given as an integer; a map is used since the indices need not
 * be contiguous.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ResponseFunctions.h,v 1.21 2003/07/21 22:14:56 jchiang Exp $
 */

class ResponseFunctions {
    
public:
    
   virtual ~ResponseFunctions() {}

   static ResponseFunctions * instance();

   void setRespPtrs(std::map<unsigned int, latResponse::Irfs *> &respPtrs)
      {m_respPtrs = respPtrs;}

   latResponse::Irfs * respPtr(unsigned int eventType);

protected:

   ResponseFunctions() {}

private:

   static ResponseFunctions * s_instance;

   std::map<unsigned int, latResponse::Irfs *> m_respPtrs;

};

} // namespace Likelihood

#endif // Likelihood_ResponseFunctions_h
