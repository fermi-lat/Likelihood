/** @file Arg.h
 * @brief Declaration of Arg class
 * $Header:
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
 * $Header: */

class Arg {
    
public:
   
   virtual ~Arg() {}

protected:

   Arg() {}

};

} // namespace Likelihood

#endif // Arg_h
