/** 
 * @file Npred.cxx
 * @brief Implementation of the Npred class, which encapsulates the
 * Npred methods of Sources in a Function context.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Npred.cxx,v 1.16 2015/03/03 18:05:37 jchiang Exp $
 */

#include <string>
#include <vector>

#include "Likelihood/Npred.h"

namespace Likelihood {

Npred::Npred() : optimizers::Function("Npred", 0, "") {
}

double Npred::value(const optimizers::Arg &x) const {
   const Source * src = dynamic_cast<const SrcArg &>(x).getValue();
   return const_cast<Source *>(src)->Npred();
}

double Npred::derivByParamImp(const optimizers::Arg & x, 
                              const std::string & paramName) const {
   const Source * src = dynamic_cast<const SrcArg &>(x).getValue();
   return const_cast<Source *>(src)->NpredDeriv(paramName);
}

void Npred::fetchDerivs(const optimizers::Arg & x,
                        std::vector<double> & derivs, 
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

void Npred::buildParameterVector(const optimizers::Arg & x) {
   m_parameter.clear();
   const Source * src = dynamic_cast<const SrcArg &>(x).getValue();

   std::vector<optimizers::Parameter> params;
   src->spectrum().getParams(params);
   for (size_t i(0); i < params.size(); i++) {
      m_parameter.push_back(params.at(i));
   }
}   

} // namespace Likelihood
