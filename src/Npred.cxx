/** 
 * @file Npred.cxx
 * @brief Implementation of the Npred class, which encapsulates the
 * Npred methods of Sources in a Function context.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/Npred.cxx,v 1.11 2010/05/18 18:27:32 jchiang Exp $
 */

#include <vector>
#include <string>
#include "Likelihood/Npred.h"

namespace Likelihood {

double Npred::value(optimizers::Arg &x) const {
   Source *src = dynamic_cast<SrcArg &>(x).getValue();

   return src->Npred();
}

double Npred::derivByParam(optimizers::Arg &x, 
                           const std::string &paramName) const {
   Source *src = dynamic_cast<SrcArg &>(x).getValue();

   double value = src->NpredDeriv(paramName);
   return value;
}

void Npred::fetchDerivs(optimizers::Arg &x, std::vector<double> &derivs, 
                        bool getFree) const {
   if (!derivs.empty()) derivs.clear();

   const_cast<Npred *>(this)->buildParameterVector(x);

   for (size_t i(0); i < m_parameter.size(); i++) {
      if (!getFree || m_parameter[i].isFree()) {
         derivs.push_back(derivByParam(x, m_parameter[i].getName()));
      }
   }
}

void Npred::buildParameterVector(optimizers::Arg & x) {
   m_parameter.clear();
   Source * src = dynamic_cast<SrcArg &>(x).getValue();

   std::vector<optimizers::Parameter> params;
   src->spectrum().getParams(params);
   for (size_t i(0); i < params.size(); i++) {
      m_parameter.push_back(params.at(i));
   }
}   

} // namespace Likelihood
