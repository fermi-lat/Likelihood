/**
 * @file MultipleScaleFactor.cxx
 * @brief User configurable multiply broken power-law.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/MultipleScaleFactor.cxx,v 1.5 2015/03/17 03:16:21 jchiang Exp $
 */

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include "optimizers/dArg.h"
#include "Likelihood/MultipleScaleFactor.h"

namespace Likelihood {
  
  MultipleScaleFactor(float e_min, float e_max, unsigned nbins)
    :optimizers::Function("MultipleScaleFactor", 100, "Normalization"),
     m_emin(e_min),
     m_emax(e_max),
     m_nbins(nbins){
    addParams();
  }
  
  MultipleScaleFactor(const MultipleScaleFactor& other)
    :optimizers::Function(other),
     m_emin(e_min),
     m_emax(e_max),
     m_nbins(nbins){
  }


  void MultipleScaleFactor::fix_params_outside_energy_range(float e_min, float e_max) {
    float log_emin = std::log10(e_min);
    float log_emax = std::log10(e_max);
    size_t i_min = find_index(log_emin);
    size_t i_max = find_index(log_emax);
    for ( size_t i(0); i < m_nbins; i++ ) {
      if ( i < i_min ) {
	fixScaleFactor(i);
      } else if ( i > i_max ) {
	fixScaleFactor(i);
      } else {
	freeScaleFactor(i);
      }
    }
  }



  size_t MultipleScaleFactor::findIndex(float logx) const {
    size_t k(std::upper_bound(m_logEnergies.begin(), m_logEnergies.end(), x)
	     - m_logEnergies.begin());
    return k;
  }

  double MultipleScaleFactor::value(const optimizers::Arg & xarg) const {

    double x(dynamic_cast<const optimizers::dArg &>(xarg).getValue());

    float logE = std::log10(x);
    size_t k = findIndex(logE);

    if ( k >= m_logEnergies.size() ) {
      throw std::runtime_error("Energy outside of range for MultipleScaleFactor");
    }

    // The vector index for the photon index parameter is offset by 1 to
    // account for the Normalization parameter.
    double value = m_parameter[0].getTrueValue();
    value * = m_parameter[k + 1].getTrueValue();
    
    return value;
  }

  double MultipleScaleFactor::derivByParamImp(const optimizers::Arg & xarg,
					      const std::string & paramName) const {

    if (paramName != "Normalization") {
      return m_parameter[0].getTrueValue();
    }

    double x(dynamic_cast<const optimizers::dArg &>(xarg).getValue());

    float logE = std::log10(x);
    size_t k = findIndex(logE);
    
    if ( k >= m_logEnergies.size() ) {
      throw std::runtime_error("Energy outside of range for MultipleScaleFactor");
    }
    
    return m_parameter[k+1].getTrueValue();
  }

  void MultipleScaleFactor::fixScaleFactor(size_t i) {
    m_parameter[i+1].setFree(false);
  }

  void MultipleScaleFactor::freeScaleFactor(size_t i) {
    m_parameter[i+1].setFree(true);
  }

  void MultipleScaleFactor::addParams() {
    addParam("Normalization", 1.0, true);
    for ( size_t i(0); i < m_nbins; i++ ) {
      std::ostringstream parname;
      parname << s_par_prefix << i;
      addParam(parname.str(), 1.0, true);
    }   
  }

} // namespace Likelihood
