/**
 * @brief Decorator class for multiplying underlying spectral models
 * by a user-selectable EBL attenuation.
 * 
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/EblAtten.cxx,v 1.1 2009/06/17 19:58:12 jchiang Exp $
 */

#include <algorithm>

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/PowerLaw2.h"
#include "Likelihood/EblAtten.h"

namespace Likelihood {

EblAtten::EblAtten() : m_spectrum(new PowerLaw2()), m_tau(0) {
   init();
}

EblAtten::EblAtten(const optimizers::Function & spectrum, 
                   double tau_norm, double redshift,
                   size_t ebl_model) 
   : m_spectrum(spectrum.clone()), m_tau(0) {
   init(tau_norm, redshift, ebl_model);
}

void EblAtten::init(double tau_norm, double redshift, size_t ebl_model) {
   setMaxNumParams(m_spectrum->getNumParams() + 3);
   
   std::vector<optimizers::Parameter> params;
   m_spectrum->getParams(params);
   for (size_t i(0); i < params.size(); i++) {
      m_parameter.push_back(params.at(i));
      m_parameter.back().setParRef(&m_spectrum->parameter(params.at(i).getName()));
   }

   addParam("tau_norm", tau_norm, true);
   addParam("redshift", redshift, false);
   addParam("ebl_model", ebl_model, false);

   m_tau = new IRB::EblAtten(static_cast<IRB::EblModel>(ebl_model));

   m_funcType = Addend;
   m_argType = "dArg";
   
   m_genericName = "EblAtten::" + m_spectrum->genericName();
   m_normParName = m_spectrum->normPar().getName();
}

EblAtten::EblAtten(const EblAtten & other) 
   : optimizers::Function(other),
     m_spectrum(other.m_spectrum->clone()),
     m_tau(new IRB::EblAtten(*other.m_tau)) {
}

EblAtten::~EblAtten() throw() {
   try {
      delete m_spectrum;
      delete m_tau;
   } catch (...) {
   }
}

double EblAtten::value(optimizers::Arg & xarg) const {
   double energy = dynamic_cast<optimizers::dArg &>(xarg).getValue();
   return m_spectrum->operator()(xarg)*attenuation(energy);
//   return m_spectrum->operator()(xarg);
}

double EblAtten::derivByParam(optimizers::Arg & xarg,
                             const std::string & paramName) const {
   double energy(dynamic_cast<optimizers::dArg &>(xarg).getValue());
   int iparam(-1);
   for (unsigned int i = 0; i < m_parameter.size(); i++) {
      if (paramName == m_parameter[i].getName()) {
         iparam = i;
	 break;
      }
   }

   if (iparam == -1) {
      throw optimizers::ParameterNotFound(paramName, getName(),
                                          m_genericName + "::derivByParam");
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
//   return m_spectrum->derivByParam(xarg, paramName);
}

// void EblAtten::setParam(const std::string & parName, double value) {
//    try {
//       m_spectrum->setParam(parName, value);
//    } catch (optimizers::ParameterNotFound &) {
//    }
//    optimizers::Function::setParam(parName, value);
// }

void EblAtten::setParam(const optimizers::Parameter & param) {
   optimizers::Parameter & my_par(parameter(param.getName()));
   my_par.setName(param.getName());
   my_par.setValue(param.getValue());
   my_par.setScale(param.getScale());
   my_par.setBounds(param.getBounds());
   my_par.setFree(param.isFree());
   my_par.setAlwaysFixed(param.alwaysFixed());
   my_par.setError(param.error());
}

// double EblAtten::getParamValue(const std::string & parName) const {
//    try {
//       return m_spectrum->getParamValue(parName);
//    } catch (optimizers::ParameterNotFound &) {
//    }
//    return optimizers::Function::getParamValue(parName);
// }

// const optimizers::Parameter & 
// EblAtten::getParam(const std::string & parName) const {
//    try {
//       return m_spectrum->getParam(parName);
//    } catch (optimizers::ParameterNotFound &) {
//    }
//    return optimizers::Function::getParam(parName);
// }

// optimizers::Parameter & EblAtten::parameter(const std::string & parName) {
//    try {
//       return m_spectrum->parameter(parName);
//    } catch (optimizers::ParameterNotFound &) {
//    }
//    return optimizers::Function::parameter(parName);
// }

// optimizers::Parameter & EblAtten::normPar() {
//    return m_spectrum->normPar();
// }

// void EblAtten::getFreeParams(std::vector<optimizers::Parameter> & pars) const {
//    m_spectrum->getFreeParams(pars);
//    for (size_t i(m_spectrum->getNumParams()); i < m_parameter.size(); i++) {
//       pars.push_back(m_parameter.at(i));
//    }
// }

// std::vector<double>::const_iterator
// EblAtten::setParamValues_(std::vector<double>::const_iterator it) {
//    it = m_spectrum->setParamValues_(it);
//    it -= m_spectrum->getNumParams();
//    return Function::setParamValues_(it);
// }

// void EblAtten::setParams(const std::vector<optimizers::Parameter> & pars) {
//    Function::setParams(pars);
//    std::vector<optimizers::Parameter> spec_pars;
//    spec_pars.resize(m_spectrum->getNumParams());
//    std::copy(pars.begin(), pars.begin() + spec_pars.size(), 
//              spec_pars.begin());
//    m_spectrum->setParams(spec_pars);
// }

// std::vector<double>::const_iterator
// EblAtten::setFreeParamValues_(std::vector<double>::const_iterator it) {
//    it = m_spectrum->setFreeParamValues_(it);
//    it -= m_spectrum->getNumFreeParams();
//    return Function::setFreeParamValues_(it);
// }

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
