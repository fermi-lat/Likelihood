#include "../Likelihood/Function.h"
#include "../Likelihood/Arg.h"

namespace Likelihood {
/** 
 * @class PowerLaw
 *
 * @brief A power-law function
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/PowerLaw.h,v 1.2 2003/02/23 22:09:48 jchiang Exp $
 */
    
class PowerLaw : public Function {
public:

   PowerLaw(){m_init(0, -2, 1);};
   PowerLaw(double Prefactor, double Index, double Scale)
      {m_init(Prefactor, Index, Scale);};

   double value(Arg&) const;

   double derivByParam(Arg &x, const std::string &paramName) const;

   double integral(Arg &xmin, Arg &xmax) const;

private:

   void m_init(double Prefactor, double Index, double Scale);

};

} // namespace Likelihood

