/** @file Rosen.h
 * @brief Declaration for a 2D Rosenbrock objective function
 * @author J. Chiang
 *
 * $Header$
 */

#include "Likelihood/Statistic.h"
#include "Likelihood/Arg.h"

namespace Likelihood {
/** 
 * @class Rosen
 *
 * @brief A 2D Rosenbrock test function
 *
 * @author J. Chiang
 *    
 * $Header$
 */
    
class Rosen : public Statistic {
public:

   Rosen() : m_prefactor(100) {init();}
   Rosen(double prefactor) : m_prefactor(prefactor) {init();}
      
   double value(const std::vector<double> &paramVec);

   void getFreeDerivs(std::vector<double> &freeDerivs);

   double derivByParam(Arg &, const std::string &paramName);

private:

   double m_prefactor;

   void init();

};

} // namespace Likelihood

