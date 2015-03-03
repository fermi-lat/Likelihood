/**
 * @brief Decorator class for restricting the underlying function
 * to a specified energy band.
 * 
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/EnergyBand.cxx,v 1.2 2011/06/10 17:18:54 jchiang Exp $
 */

#include <algorithm>
#include <iostream>

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/PowerLaw2.h"
#include "Likelihood/EnergyBand.h"

namespace Likelihood {

EnergyBand::EnergyBand() 
   : optimizers::Function("EnergyBand::PowerLaw2", 6, "Integral"),
     m_spectrum(new PowerLaw2()) {
   init(100, 3e5);
}

EnergyBand::EnergyBand(const optimizers::Function & spectrum, 
                       double emin, double emax)
   : optimizers::Function("EnergyBand::" + spectrum.genericName(),
                          spectrum.getNumParams() + 2,
                          spectrum.normPar().getName()),
     m_spectrum(spectrum.clone()) {
   init(emin, emax);
}
   
void EnergyBand::init(double emin, double emax) {
   std::vector<optimizers::Parameter> params;
   m_spectrum->getParams(params);
   for (size_t i(0); i < params.size(); i++) {
      m_parameter.push_back(params[i]);
   }
   setParRefs();

   bool allowed_free;
   addParam("Emin", emin, allowed_free=false);
   addParam("Emax", emax, allowed_free=false);
}

void EnergyBand::setParRefs() {
   for (size_t i(0); i < m_spectrum->getNumParams(); i++) {
      std::string name(m_parameter[i].getName());
      m_parameter[i].setParRef(&m_spectrum->parameter(name));
   }
}

EnergyBand::EnergyBand(const EnergyBand & other) 
   : optimizers::Function(other),
     m_spectrum(other.m_spectrum->clone()) {
   setParRefs();
}

EnergyBand & EnergyBand::operator=(const EnergyBand & rhs) {
   if (this != &rhs) {
      m_spectrum = rhs.m_spectrum->clone();
      setParRefs();
   }
   return *this;
}

EnergyBand::~EnergyBand() throw() {
   try {
      delete m_spectrum;
   } catch (...) {
   }
}
   
double EnergyBand::value(optimizers::Arg & xarg) const {
   double energy = dynamic_cast<optimizers::dArg &>(xarg).getValue();
   int emin(m_spectrum->getNumParams());
   int emax(m_spectrum->getNumParams() + 1);
   
   if (energy < m_parameter[emin].getTrueValue() ||
       energy >= m_parameter[emax].getTrueValue()) {
      return 0;
   }
   return m_spectrum->operator()(xarg);
}
   
double EnergyBand::derivByParamImp(optimizers::Arg & xarg,
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
                                          genericName() + "::derivByParam");
   }
   
   int emin(m_spectrum->getNumParams());
   int emax(m_spectrum->getNumParams() + 1);

   if (iparam >= emin) {
      throw std::runtime_error("EnergyBand: Attempt to take derivative wrt "
                               "a fixed parameter.");
   }
   
   if (energy >= m_parameter[emin].getTrueValue() &&
       energy < m_parameter[emax].getTrueValue()) {
      return m_spectrum->derivByParam(xarg, paramName);
   } 

   return 0;
}

void EnergyBand::setParam(const optimizers::Parameter & param) {
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
