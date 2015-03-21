/**
 * @brief Decorator class for multiplying underlying spectral models
 * by a user-selectable EBL attenuation.
 * 
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/EblAtten.cxx,v 1.6 2015/03/03 18:05:37 jchiang Exp $
 */

#include <algorithm>
#include <iostream>

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/PowerLaw2.h"
#include "Likelihood/EblAtten.h"

namespace Likelihood {

EblAtten::EblAtten() 
   : optimizers::Function("EblAtten::PowerLaw2", 7, "Integral"),
     m_spectrum(new PowerLaw2()), m_tau(0) {
   init();
}

EblAtten::EblAtten(const optimizers::Function & spectrum, 
                   double tau_norm, double redshift,
                   size_t ebl_model) 
   : optimizers::Function("EblAtten::" + spectrum.genericName(),
                          spectrum.getNumParams() + 3, 
                          spectrum.normPar().getName()),
     m_spectrum(spectrum.clone()), m_tau(0) {
   init(tau_norm, redshift, ebl_model);
}

void EblAtten::init(double tau_norm, double redshift, size_t ebl_model) {
   std::vector<optimizers::Parameter> params;
   m_spectrum->getParams(params);
   for (size_t i(0); i < params.size(); i++) {
      m_parameter.push_back(params.at(i));
   }
   setParRefs();

   addParam("tau_norm", tau_norm, true);
   addParam("redshift", redshift, false);
   addParam("ebl_model", ebl_model, false);

   m_tau = new IRB::EblAtten(static_cast<IRB::EblModel>(ebl_model));
}

void EblAtten::setParRefs() {
   for (size_t i(0); i < m_spectrum->getNumParams(); i++) {
      std::string name(m_parameter.at(i).getName());
      m_parameter.at(i).setParRef(&m_spectrum->parameter(name));
   }
}

EblAtten::EblAtten(const EblAtten & other) 
   : optimizers::Function(other),
     m_spectrum(other.m_spectrum->clone()),
     m_tau(new IRB::EblAtten(*other.m_tau)) {
   setParRefs();
}

EblAtten & EblAtten::operator=(const EblAtten & rhs) {
   if (this != &rhs) {
      delete m_spectrum;
      delete m_tau;
      m_spectrum = rhs.m_spectrum->clone();
      m_tau = new IRB::EblAtten(*rhs.m_tau);
      setParRefs();
   }
   return *this;
}

EblAtten::~EblAtten() throw() {
   try {
      delete m_spectrum;
      delete m_tau;
   } catch (...) {
   }
}

double EblAtten::value(const optimizers::Arg & xarg) const {
   double energy = dynamic_cast<const optimizers::dArg &>(xarg).getValue();
   return m_spectrum->operator()(xarg)*attenuation(energy);
}

double EblAtten::derivByParamImp(const optimizers::Arg & xarg,
                                 const std::string & paramName) const {
   double energy(dynamic_cast<const optimizers::dArg &>(xarg).getValue());
   int iparam(-1);
   for (unsigned int i = 0; i < m_parameter.size(); i++) {
      if (paramName == m_parameter[i].getName()) {
         iparam = i;
	 break;
      }
   }

   if (iparam == -1) {
      throw optimizers::ParameterNotFound(paramName, getName(),
                                          genericName() + "::derivByParam");
   }

   if (iparam < static_cast<int>(m_spectrum->getNumParams())) {
      return attenuation(energy)*m_spectrum->derivByParam(xarg, paramName);
   }

   int tau_norm(m_spectrum->getNumParams());
   int redshift(tau_norm + 1);

   if (iparam != tau_norm) {
      throw std::runtime_error("EblAtten: Attempt to take derivative wrt "
                               "a fixed parameter.");
   }
   double zz(m_parameter[redshift].getTrueValue());
   return -value(xarg)*m_tau->operator()(energy, zz);
}

void EblAtten::setParam(const optimizers::Parameter & param) {
// This version preserves the references to the m_spectrum parameters.
   optimizers::Parameter & my_par(parameter(param.getName()));
   my_par.setName(param.getName());
   my_par.setValue(param.getValue());
   my_par.setScale(param.getScale());
   my_par.setBounds(param.getBounds());
   my_par.setFree(param.isFree());
   my_par.setAlwaysFixed(param.alwaysFixed());
   my_par.setError(param.error());
}

double EblAtten::attenuation(double energy) const {
   int tau_norm(m_spectrum->getNumParams());
   int redshift(tau_norm + 1);
   int ebl_model(redshift + 1);
   
   IRB::EblModel new_model = static_cast<IRB::EblModel>(static_cast<int>(m_parameter[ebl_model].getTrueValue()));
   if (new_model != m_tau->model()) {
      delete m_tau;
      m_tau = new IRB::EblAtten(new_model);
   }

   double zz(m_parameter[redshift].getTrueValue());
   double atten(std::exp(-m_parameter[tau_norm].getTrueValue()
                         *m_tau->operator()(energy, zz)));
   return atten;
}

} // namespace Likelihood
