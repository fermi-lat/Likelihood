/** 
 * @file PowerLaw2.h
 * @brief Declaration for the PowerLaw2 Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/PowerLaw2.h,v 1.5 2015/03/21 05:38:03 jchiang Exp $
 */

#ifndef Likelihood_PowerLaw2_h
#define Likelihood_PowerLaw2_h

#include "optimizers/Arg.h"
#include "optimizers/Function.h"

namespace Likelihood {

/** 
 * @class PowerLaw2
 *
 * @brief A power-law function that uses integrated flux and index
 * as free parameters and upper and lower bounds of integration as 
 * fixed parameters.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/PowerLaw2.h,v 1.5 2015/03/21 05:38:03 jchiang Exp $
 */
    
class PowerLaw2 : public optimizers::Function {

public:

   PowerLaw2(double Integral=1., double Index=-2., 
             double LowerLimit=100., double UpperLimit=2e5);

   double integral(optimizers::Arg & xmin, optimizers::Arg & xmax) const;

   virtual Function * clone() const {
      return new PowerLaw2(*this);
   }

protected:
   double value(const optimizers::Arg & x) const;

   double derivByParamImp(const optimizers::Arg & x, 
                          const std::string & paramName) const;

private:

   void updateCache(const double x,
		    const double gamma, const double xlo, const double xhi) 
     const {
     if(gamma != m_cGamma || xlo != m_cXLo || xhi != m_cXHi)
       {
	 double one_p_gamma = 1.+gamma;
	 m_cGamma  = gamma;
	 m_cXLo    = xlo;
	 m_cXHi    = xhi;
	 m_cLogXLo = std::log(xlo);
	 m_cLogXHi = std::log(xhi);
	 m_cPowXLo = std::pow(xlo,one_p_gamma);
	 m_cPowXHi = std::pow(xhi,one_p_gamma);
         if (gamma != -1.0) {
            m_cGXFact = one_p_gamma/(m_cPowXHi - m_cPowXLo);
         }
	 m_cX      = x + 1.0; // force recalculation of m_cPowX below
       } 

     if(s_cX != x)
       {
	 s_cX      = x;
	 s_cLogX   = std::log(x);
       }

     if(m_cX != x)
       {
	 m_cX      = x;
	 m_cPowX   = std::exp(s_cLogX*gamma);
       }
   }

   // cached variables
   mutable double m_cGamma;
   mutable double m_cXLo;
   mutable double m_cXHi;
   mutable double m_cLogXLo;
   mutable double m_cLogXHi;
   mutable double m_cPowXLo;
   mutable double m_cPowXHi;
   mutable double m_cGXFact;
   mutable double m_cX;
   mutable double m_cPowX;
   static  double s_cX;
   static  double s_cLogX;
  };

} // namespace Likelihood

#endif // Likelihood_PowerLaw2_h
