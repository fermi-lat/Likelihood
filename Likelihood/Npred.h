/** 
 * @file Npred.h
 * @brief Declaration of Npred class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Npred.h,v 1.5 2003/06/11 17:08:02 jchiang Exp $
 */

#ifndef Likelihood_Npred_h
#define Likelihood_Npred_h

#include "Likelihood/Function.h"
#include "Likelihood/Source.h"
#include "Likelihood/SrcArg.h"

namespace Likelihood {

/** 
 * @class Npred
 *
 * @brief This class encapsulates the Npred methods of Sources in a
 * Function context.
 *  
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Npred.h,v 1.5 2003/06/11 17:08:02 jchiang Exp $
 */

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

#endif // Likelihood_Npred_h
