#include <vector>
#include <string>
#include "../Likelihood/Npred.h"

namespace Likelihood {

double Npred::value(Arg &x) const {
   Source *src = dynamic_cast<SrcArg &>(x).getValue();

   return src->Npred();
}

double Npred::derivByParam(Arg &x, const std::string &paramName) const {
   Source *src = dynamic_cast<SrcArg &>(x).getValue();

   double value = src->NpredDeriv(paramName);
   return value;
}

void Npred::fetchDerivs(Arg &x, std::vector<double> &derivs, 
                        bool getFree) const {
   if (!derivs.empty()) derivs.clear();

   buildParameterVector(x);

   for (unsigned int i = 0; i < m_parameter.size(); i++) {
      if (!getFree || m_parameter[i].isFree())
         derivs.push_back(derivByParam(x, m_parameter[i].getName()));
   }
}

void Npred::buildParameterVector(Arg &x) const {
   m_parameter.clear();
   Source *src = dynamic_cast<SrcArg &>(x).getValue();

   Source::FuncMap srcFuncs = src->getSrcFuncs();
   Source::FuncMap::iterator func_it = srcFuncs.begin();
   for (; func_it != srcFuncs.end(); func_it++) {
      std::vector<Parameter> params;
      (*func_it).second->getParams(params);
      for (unsigned int i = 0; i < params.size(); i++)
         m_parameter.push_back(params[i]);
   }
}   


} // namespace Likelihood
