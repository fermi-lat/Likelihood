/** 
 * @file Npred.h
 * @brief Declaration of Npred class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Npred.h,v 1.9 2004/12/22 06:06:47 jchiang Exp $
 */

#ifndef Likelihood_Npred_h
#define Likelihood_Npred_h

#include "optimizers/Function.h"
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Npred.h,v 1.9 2004/12/22 06:06:47 jchiang Exp $
 */

class Npred : public optimizers::Function {
    
public:

   Npred() {m_genericName = "Npred";}
   virtual ~Npred() {}

   double value(optimizers::Arg &) const;
   double derivByParam(optimizers::Arg &, const std::string &) const;

protected:

   Npred * clone() const {
      return new Npred(*this);
   }

private:

   void fetchDerivs(optimizers::Arg &, std::vector<double> &derivs, 
                    bool getFree) const;
   void buildParameterVector(optimizers::Arg &);

};

} // namespace Likelihood

#endif // Likelihood_Npred_h
