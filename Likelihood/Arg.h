/** 
 * @file Arg.h
 * @brief Declaration of Arg class
 * @author J. Chiang
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Arg.h,v 1.5 2003/05/29 00:29:39 jchiang Exp $
 */

#ifndef Likelihood_Arg_h
#define Likelihood_Arg_h

namespace Likelihood {

/** 
 * @class Arg
 *
 * @brief An abstract class that encapsulates argument type
 * information so that Function's value() and Parameter passing
 * methods can be overloaded transparently.
 *
 * @authors J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Arg.h,v 1.5 2003/05/29 00:29:39 jchiang Exp $
 */

class Arg {
    
public:

   Arg() {}
   
   virtual ~Arg() {}

protected:

//   Arg() {}

};

} // namespace Likelihood

#endif // Likelihood_Arg_h
