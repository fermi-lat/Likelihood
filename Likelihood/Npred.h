/** @file Npred.h
 * @brief Declaration of Npred class
 * $Header:
 */

#ifndef Npred_h
#define Npred_h

#include "../Likelihood/Function.h"
#include "../Likelihood/Source.h"
#include "../Likelihood/Arg.h"

namespace Likelihood {

/** 
 * @class Npred
 *
 * @brief This class encapsulates the Npred methods of Sources in a
 * Function context.
 *
 * Here Arg is used only as place-holder.
 *  
 * @author J. Chiang
 *    
 * $Header: */

class Npred : public Function {
    
public:

   Npred(Source *src) : m_source(src) {}
   virtual ~Npred() {}

   double value(Arg &) const;
   double derivByParam(Arg &, const std::string &) const;

private:

   Source *m_source;

   void fetchDerivs(Arg &, std::vector<double> &derivs, bool getFree) const;

};

} // namespace Likelihood

#endif // Npred_h
