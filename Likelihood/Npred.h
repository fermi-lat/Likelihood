/** 
 * @file Npred.h
 * @brief Declaration of Npred class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/Npred.h,v 1.10 2005/02/15 00:34:42 jchiang Exp $
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
 */

class Npred : public optimizers::Function {
    
public:

   Npred() : m_use_ebounds(false), m_emin(0), m_emax(0) {
      m_genericName = "Npred";
   }

   virtual ~Npred() {}

   double value(optimizers::Arg &) const;

   double derivByParam(optimizers::Arg &, const std::string &) const;

   void set_ebounds(double emin, double emax);

   void unset_ebounds();

protected:

   Npred * clone() const {
      return new Npred(*this);
   }

private:

   bool m_use_ebounds;
   double m_emin;
   double m_emax;

   void fetchDerivs(optimizers::Arg &, std::vector<double> &derivs, 
                    bool getFree) const;
   void buildParameterVector(optimizers::Arg &);

};

} // namespace Likelihood

#endif // Likelihood_Npred_h
