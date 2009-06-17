/**
 * @brief Decorator class for multiplying underlying spectral models
 * by a user-selectable EBL attenuation.
 * 
 * @author J. Chiang
 * 
 * $Header$
 */

//namespace Likelihood {

EblAtten::EblAtten(const optimizers::Function & spectrum, 
                   double tau_norm, double redshift, IRB::EblModel model) 
   : m_spectrum(spectrum.clone()), m_tau(0) {
   setMaxNumParams(m_spectrum->getNumParams() + 3);
   
   std::vector<optimizers::Parameter> params;
   m_spectrum->getParams(params);
   for (size_t i(0); i < params.size(); i++) {
      addParam(params.at(i));
   }

   addParam("tau_norm", redshift, true);
   addParam("redshift", redshift, false);
   addParam("ebl_model", ebl_model, false);

   m_tau = new IRB::EblAtten(model);

   m_funcType = Addend;
   m_argType = "dArg";
   
   m_genericName = "EblAtten::" + m_spectrum->genericName();
   m_normParName = m_spectrum->normPar().getName();
}

double EblAtten::value(optimizers::Arg & xarg) const {
   double energy = dynamic_cast<optimizers::dArg &>(xarg).getValue();
   return m_spectrum->operator()(xarg)*attenuation(energy);
}

double EbAtten::derivByParam(optimizers::Arg & xarg,
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

   if (iparam < m_spectrum->getNumParams()) {
      return attenuation(energy)*m_spectrum->derivByParam(xarg, paramName);
   }

   enum ParamTypes {tau_norm=m_spectrum->getNumParams(), redshift, ebl_model};

   if (iparam != tau_norm) {
      throw std::runtime_error("EblAtten: Attempt to take derivative wrt "
                               "a fixed parameter.");
   }
   double zz(m_parameter[redshift].getTrueValue());
   return value(xarg)*m_tau->operator()(energy, zz));
}

void EblAtten::setParam(const std::string & parName, double value) {
   try {
      m_spectrum->setParam(parName, value);
   } catch (optimizers::ParameterNotFound &) {
      optimizers::Function::setParam(parName, value);
   }
}

void EblAtten::setParam(const optimizers::Parameter & param) {
   try {
      m_spectrum->setParam(parName, value);
   } catch (optimizers::ParameterNotFound &) {
      optimizers::Function::setParam(parName, value);
   }
}

double EblAtten::getParamValue(const std::string & parName) const {
   try {
      return m_spectrum->getParamValue(parName);
   } catch (optimizers::ParameterNotFound &) {
   }
   return optimizers::Function::getParamValue(parName);
}

optimizers::Parameter & EblAtten::getParam(const std::string & parName) const {
   try {
      return m_spectrum->getParam(parName);
   } catch (optimizers::ParameterNotFound &) {
   }
   return optimizers::Function::getParam(parName);
}

optimizers::Parameter & EblAtten::normPar() {
   return m_spectrum->normPar();
}

void 


double EblAtten::attenuation(double energy) const {
   enum ParamTypes {tau_norm=m_spectrum->getNumParams(), redshift, ebl_model};

   if (m_parameter[ebl_model] != m_tau->model()) {
      delete m_tau;
      m_tau = new IRB::EblAtten(m_parameter[ebl_model]);
   }

   double zz(m_parameter[redshift].getTrueValue());
   double atten(std::exp(-m_parameter["tau_norm"].getTrueValue()
                         *m_tau->operator()(energy, zz)));
   return atten;
}

bool is_ebl_param(const std::string & paramName) const {
   try {
      m_spectrum->getParam(paramName);
   } catch (optimizers::ParameterNotFound &) {
      return true;
   }
   return false;
}

} // namespace Likelihood
