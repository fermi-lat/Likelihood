#include <vector>
#include <string>
#include "../Likelihood/Npred.h"

namespace Likelihood {

double Npred::value(Arg &x) const {
// Note that Arg is a dummy variable that is not used by Npred
// directly except as a place-holder for Function's derivative passing
// methods.

   return m_source->Npred();
}

double Npred::derivByParam(Arg &x, const std::string &paramName) const {

   return m_source->NpredDeriv(paramName);
}

void Npred::fetchDerivs(Arg &x, std::vector<double> &derivs, 
                        bool getFree) const {
   if (!derivs.empty()) derivs.clear();

   Source::FuncMap srcFuncs = m_source->getSrcFuncs();
   Source::FuncMap::iterator func_it = srcFuncs.begin();
   for (; func_it != srcFuncs.end(); func_it++) {
      std::vector<double> my_derivs;
      if (getFree) {
         (*func_it).second->getFreeDerivs(x, my_derivs);
      } else {
         (*func_it).second->getDerivs(x, my_derivs);
      }
      for (unsigned int i = 0; i < my_derivs.size(); i++) 
         derivs.push_back(my_derivs[i]);
   }
}

} // namespace Likelihood
