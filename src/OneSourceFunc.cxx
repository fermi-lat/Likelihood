/**
 * @file OneSourceFunc.cxx
 * @brief Extended likelihood function for one source.
 *
 * @author P. Nolan
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/OneSourceFunc.cxx,v 1.5 2005/03/04 22:08:27 jchiang Exp $
 */

#include "Likelihood/OneSourceFunc.h"
#include "Likelihood/Source.h"
#include "optimizers/Parameter.h"
#include <cmath>

namespace Likelihood {

  //! cutoff values for minimum weight and flux

  OneSourceFunc::OneSourceFunc(Source * src,
			       const std::vector<Event>& evt,
			       std::vector<double> * weights):
    Statistic("OneSourceFunc", 0),
    m_src(src),
    m_events(evt),
    m_weights(weights),
    m_epsw(1.e-3),
    m_epsf(1.e-20)
  {
     setName("OneSourceFunc");
    syncParams();
  }

  void OneSourceFunc::setEpsW(double e) {
    m_epsw = e;
  }

  void OneSourceFunc::setEpsF(double e) {
    m_epsf = e;
  }

  double OneSourceFunc::value(optimizers::Arg& ) const {

    double val = 0.;
    //    double wtot = 0.;
    //    int nused = 0;
    for (unsigned int i = 0; i < m_events.size(); i++) {
      double w = (m_weights == NULL) ? 1.0 : (*m_weights)[i];
      if (w > m_epsw) {
	double q = m_src->fluxDensity(m_events[i]);
	if (fabs(q) > m_epsf) {
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

  double OneSourceFunc::derivByParamImp(optimizers::Arg &, 
                                        const std::string & paramName) const {

    double deriv = 0;
    //    double wtot = 0.;
    for (unsigned int i = 0; i < m_events.size(); i++) {
      double w = (m_weights == NULL) ? 1.0 : (*m_weights)[i];
      if (w > m_epsw) {
	//	wtot += w;
	double q = m_src->fluxDensity(m_events[i]);
	if (fabs(q) > m_epsf) {
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
