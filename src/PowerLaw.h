#include "../Likelihood/Function.h"

namespace Likelihood {
/** 
 * @class PowerLaw
 *
 * @brief A power-law function
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools/Likelihood/src/Function.h,v 1.1.1.1 2003/01/30 23:23:03 burnett Exp $
 */
    
class PowerLaw : public Function {
public:

   PowerLaw(){m_init(0, -2, 1);};
   PowerLaw(const double Prefactor, const double Index, const double Scale)
      {m_init(Prefactor, Index, Scale);};
   PowerLaw(const PowerLaw&);

   virtual double value(const double) const;
   virtual double operator()(const double x) const {return value(x);};
   virtual double derivByParam(const double, 
			       const std::string paramName) const;
   virtual std::vector<double> getDerivs(const double) const;
   virtual std::vector<double> getFreeDerivs(const double) const;

private:

   void m_init(const double Prefactor, const double Index, const double Scale);

};

} // namespace Likelihood

