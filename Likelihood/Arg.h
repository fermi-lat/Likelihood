/** 
 * @file Arg.h
 * @brief Declaration of Arg class
 * @author J. Chiang
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Arg.h,v 1.4 2003/03/17 00:53:42 jchiang Exp $
 */

#ifndef Arg_h
#define Arg_h

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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Arg.h,v 1.4 2003/03/17 00:53:42 jchiang Exp $
 */

class Arg {
    
public:

   Arg() {}
   
   virtual ~Arg() {}

protected:

//   Arg() {}

};

} // namespace Likelihood

#endif // Arg_h
