/**
 * @file MultipleScaleFactor.h
 * @brief User configurable set of energy dependent scale factors
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/MultipleBrokenPowerLaw.h,v 1.2 2015/03/03 18:05:36 jchiang Exp $
 */

#ifndef Likelihood_MultpleScaleFactor_h
#define Likelihood_MultpleScaleFactor_h

#include "optimizers/Arg.h"
#include "optimizers/Function.h"

namespace Likelihood {

  /**
   * @class MultipleScaleFactor
   */
  
  class MultipleScaleFactor : public optimizers::Function {
    
  private:

    static const std::string s_par_prefix;

  public:
    
    MultipleScaleFactor(float e_min, float e_max, unsigned nbins);
  
    MultipleScaleFactor(const MultipleScaleFactor& other);

    virtual Function * clone() const {
      return new MultipleScaleFactor(*this);
    }

    void fix_params_outside_energy_range(float e_min, float e_max);
    
  protected:
    
    size_t findIndex(float logx) const;

    double value(const optimizers::Arg & x) const;
    
    double derivByParamImp(const optimizers::Arg & x,
			   const std::string & paramName) const;

    void fixScaleFactor(size_t i);

    void freeScaleFactor(size_t i);

  private:

    void addParams();
    
    const float m_emin;
    const float m_emax;
    const size_t m_nbins;

    std::vector<float> m_logEnergies;
    
  };
  
} // namespace Likelihood 

#endif // Likelihood_MultpleScaleFactor_h
