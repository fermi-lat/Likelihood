#include "Likelihood/OneSourceFunc.h"
#include "Likelihood/Source.h"
#include "optimizers/Parameter.h"
#include <cmath>

namespace Likelihood {

  //! cutoff values for minimum weight and flux
  const double epsw = 1.e-3;
  const double epsf = 1.e-20;

  OneSourceFunc::OneSourceFunc(Source * src,
			       std::vector<Event>& evt,
			       std::vector<double> * weights):
    m_src(src),
    m_events(evt),
    m_weights(weights) 
  {
    m_functionName = "OneSourceFunc";
    syncParams();
  }

  double OneSourceFunc::value(optimizers::Arg& arg) const {

    double val = 0.;
    //    double wtot = 0.;
    //    int nused = 0;
    for (unsigned int i = 0; i < m_events.size(); i++) {
      double w = (m_weights == NULL) ? 1.0 : (*m_weights)[i];
      if (w > epsw) {
	double q = m_src->fluxDensity(m_events[i]);
	if (abs(q) > epsf) {
	  val += w * log(q);
	  //	  wtot += w;
	  //	  nused++;
	}
      }
    }
    double bar = m_src->Npred();
    double foo = val - bar;
    //  std::cout << "OSF:v value " << foo << " photons " << nused << " wtot " 
    //	      << wtot <<  " sum " << val
    //	      << " integral " << bar << std::endl;
    return foo;
  }

  double OneSourceFunc::derivByParam(optimizers::Arg& x, 
				     const std::string& paramName) const {

    double deriv = 0;
    //    double wtot = 0.;
    for (unsigned int i = 0; i < m_events.size(); i++) {
      double w = (m_weights == NULL) ? 1.0 : (*m_weights)[i];
      if (w > epsw) {
	//	wtot += w;
	double q = m_src->fluxDensity(m_events[i]);
	if (abs(q) > epsf) {
	  double v = m_src->fluxDensityDeriv(m_events[i], paramName);
	  deriv += w * v / q;
	}
      }
    }
    double bar = m_src->NpredDeriv(paramName);
    double foo = deriv - bar;
    return foo;
  }

  void OneSourceFunc::syncParams(void) {
    m_parameter.clear();
    Source::FuncMap srcFuncs = m_src->getSrcFuncs();
    Source::FuncMap::iterator func_it = srcFuncs.begin();
    for (; func_it != srcFuncs.end(); func_it++) {
      std::vector<optimizers::Parameter> params;
      (*func_it).second->getParams(params);
      // m_parameter.insert(m_parameter.end(), params.begin(), params.end());
      for (unsigned int ip = 0; ip < params.size(); ip++ )
	m_parameter.push_back(params[ip]);
    }
  }

  std::vector<double>::const_iterator 
  OneSourceFunc::setFreeParamValues_(std::vector<double>::const_iterator it) {
    syncParams();
    Likelihood::Source::FuncMap srcFuncs = m_src->getSrcFuncs();
    Likelihood::Source::FuncMap::iterator func_it = srcFuncs.begin();
    for (; func_it != srcFuncs.end(); func_it++) {
      it = (*func_it).second->setFreeParamValues_(it);
    }
    syncParams();
    return it;
  }

  std::vector<double>::const_iterator 
  OneSourceFunc::setParamValues_(std::vector<double>::const_iterator it) {
    syncParams();
    Likelihood::Source::FuncMap srcFuncs = m_src->getSrcFuncs();
    Likelihood::Source::FuncMap::iterator func_it = srcFuncs.begin();
    for (; func_it != srcFuncs.end(); func_it++) {
      it = (*func_it).second->setParamValues_(it);
    }
    syncParams();
    return it;
  }

  void OneSourceFunc::setParams(std::vector<optimizers::Parameter> &params)
    throw (optimizers::Exception, optimizers::ParameterNotFound) {
    unsigned int numParams = getNumParams();
    if (params.size() != numParams) 
      throw Exception
	("OneSourceFunc::setParams: Inconsistent number of Parameters");
    int k = 0;
    Source::FuncMap srcFuncs = m_src->getSrcFuncs();
    Source::FuncMap::iterator func_it = srcFuncs.begin();
    for (; func_it != srcFuncs.end(); func_it++) {
      unsigned int numParams = func_it->second->getNumParams();
      for (unsigned int j = 0; j < numParams; j++, k++)
	func_it->second->setParam(params[k]);
    }
    syncParams();
  }

}
