/** 
 * @file Npred.cxx
 * @brief Implementation of the Npred class, which encapsulates the
 * Npred methods of Sources in a Function context.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Npred.cxx,v 1.8 2004/09/15 23:12:37 jchiang Exp $
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

   for (unsigned int i = 0; i < m_parameter.size(); i++) {
      if (!getFree || m_parameter[i].isFree())
         derivs.push_back(derivByParam(x, m_parameter[i].getName()));
   }
}

//void Npred::buildParameterVector(optimizers::Arg &x) const {
void Npred::buildParameterVector(optimizers::Arg &x) {
   m_parameter.clear();
   Source *src = dynamic_cast<SrcArg &>(x).getValue();

   Source::FuncMap srcFuncs = src->getSrcFuncs();
   Source::FuncMap::const_iterator func_it = srcFuncs.begin();
   for (; func_it != srcFuncs.end(); ++func_it) {
      std::vector<optimizers::Parameter> params;
      (*func_it).second->getParams(params);
      for (unsigned int i = 0; i < params.size(); i++)
         m_parameter.push_back(params[i]);
   }
}   

} // namespace Likelihood
