/** @file Npred.h
 * @brief Declaration of Npred class
 * $Header:
 */

#ifndef Npred_h
#define Npred_h

#include "../Likelihood/Function.h"
#include "../Likelihood/Source.h"
#include "../Likelihood/SrcArg.h"

namespace Likelihood {

/** 
 * @class Npred
 *
 * @brief This class encapsulates the Npred methods of Sources in a
 * Function context.
 *  
 * @author J. Chiang
 *    
 * $Header: */

class Npred : public Function {
    
public:

   Npred() {}
   virtual ~Npred() {}

   double value(Arg &) const;
   double derivByParam(Arg &, const std::string &) const;

private:

   void fetchDerivs(Arg &, std::vector<double> &derivs, bool getFree) const;
   void buildParameterVector(Arg &) const;

};

} // namespace Likelihood

#endif // Npred_h
