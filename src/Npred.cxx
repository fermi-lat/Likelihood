/** 
 * @file Npred.cxx
 * @brief Implementation of the Npred class, which encapsulates the
 * Npred methods of Sources in a Function context.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Npred.cxx,v 1.15 2012/06/27 20:31:51 jchiang Exp $
 */

#include <string>
#include <vector>

#include "Likelihood/Npred.h"

namespace Likelihood {

Npred::Npred() : optimizers::Function("Npred", 0, "") {
}

double Npred::value(optimizers::Arg &x) const {
   Source * src = dynamic_cast<SrcArg &>(x).getValue();
   return src->Npred();
}

double Npred::derivByParamImp(optimizers::Arg & x, 
                              const std::string & paramName) const {
   Source * src = dynamic_cast<SrcArg &>(x).getValue();
   return src->NpredDeriv(paramName);
}

void Npred::fetchDerivs(optimizers::Arg & x, std::vector<double> & derivs, 
                        bool getFree) const {
   if (!derivs.empty()) {
      derivs.clear();
   }

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
