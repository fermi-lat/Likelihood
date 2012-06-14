/** 
 * @file Npred.cxx
 * @brief Implementation of the Npred class, which encapsulates the
 * Npred methods of Sources in a Function context.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/Npred.cxx,v 1.13 2012/06/14 02:01:25 jchiang Exp $
 */

#include <string>
#include <vector>

#include "Likelihood/Npred.h"

namespace Likelihood {

double Npred::value(optimizers::Arg &x) const {
   Source * src = dynamic_cast<SrcArg &>(x).getValue();

   if (m_use_ebounds) {
      return src->Npred(m_emin, m_emax);
   }
   return src->Npred();
}

double Npred::derivByParam(optimizers::Arg & x, 
                           const std::string & paramName) const {
   Source * src = dynamic_cast<SrcArg &>(x).getValue();
   double value(0);
   if (m_use_ebounds) {
      value = src->NpredDeriv(paramName, m_emin, m_emax);
   } else {
      value = src->NpredDeriv(paramName);
   }
   return value;
}

void Npred::set_ebounds(double emin, double emax) {
   m_use_ebounds = true;
   m_emin = emin;
   m_emax = emax;
}

void Npred::unset_ebounds() {
   m_use_ebounds = false;
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
