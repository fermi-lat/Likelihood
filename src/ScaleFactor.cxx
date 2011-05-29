/**
 * @brief Decorator class for multiplying underlying spectral models
 * by an overall scale factor (on top of the usual normalization
 * parameter.)  Requested by DM group.
 * 
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/ScaleFactor.cxx,v 1.4 2009/06/19 06:36:14 jchiang Exp $
 */

#include <algorithm>
#include <iostream>

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/PowerLaw2.h"
#include "Likelihood/ScaleFactor.h"

namespace Likelihood {

ScaleFactor::ScaleFactor() : m_spectrum(new PowerLaw2()) {
   init();
}

ScaleFactor::ScaleFactor(const optimizers::Function & spectrum, 
                         double scale_factor)
   : m_spectrum(spectrum.clone()) {
   init(scale_factor);
}
   
void ScaleFactor::init(double scale_factor) {
   setMaxNumParams(m_spectrum->getNumParams() + 1);
   
   std::vector<optimizers::Parameter> params;
   m_spectrum->getParams(params);
   for (size_t i(0); i < params.size(); i++) {
      m_parameter.push_back(params[i]);
   }
   setParRefs();

   addParam("ScaleFactor", scale_factor, true);

   m_funcType = Addend;
   m_argType = "dArg";
   
   m_genericName = "ScaleFactor::" + m_spectrum->genericName();
   m_normParName = m_spectrum->normPar().getName();
}

void ScaleFactor::setParRefs() {
   for (size_t i(0); i < m_spectrum->getNumParams(); i++) {
      std::string name(m_parameter[i].getName());
      m_parameter[i].setParRef(&m_spectrum->parameter(name));
   }
}

ScaleFactor::ScaleFactor(const ScaleFactor & other) 
   : optimizers::Function(other),
     m_spectrum(other.m_spectrum->clone()) {
   setParRefs();
}

ScaleFactor & ScaleFactor::operator=(const ScaleFactor & rhs) {
   if (this != &rhs) {
      m_spectrum = rhs.m_spectrum->clone();
      setParRefs();
   }
   return *this;
}

ScaleFactor::~ScaleFactor() throw() {
   try {
      delete m_spectrum;
   } catch (...) {
   }
}
   
double ScaleFactor::value(optimizers::Arg & xarg) const {
   double energy = dynamic_cast<optimizers::dArg &>(xarg).getValue();
   int scale_factor(m_spectrum->getNumParams());
   return m_spectrum->operator()(xarg)*m_parameter[scale_factor].getTrueValue();
}
   
double ScaleFactor::derivByParam(optimizers::Arg & xarg,
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
   
   int scale_factor(m_spectrum->getNumParams());
   if (iparam < scale_factor) {
      return (m_parameter[scale_factor].getTrueValue()*
              m_spectrum->derivByParam(xarg, paramName));
   }

   return m_spectrum->operator()(xarg)*m_parameter[scale_factor].getScale();
}

void ScaleFactor::setParam(const optimizers::Parameter & param) {
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

} // namespace Likelihood
