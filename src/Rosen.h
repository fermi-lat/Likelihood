/** 
 * @file Rosen.h
 * @brief Declaration for a 2D Rosenbrock objective function
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Rosen.h,v 1.4 2003/06/03 23:54:23 jchiang Exp $
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Rosen.h,v 1.4 2003/06/03 23:54:23 jchiang Exp $
 */
    
class Rosen : public Statistic {
public:

   Rosen() : m_prefactor(100) {init();}
   Rosen(double prefactor) : m_prefactor(prefactor) {init();}
      
   double value(const std::vector<double> &paramVec);

   void getFreeDerivs(std::vector<double> &freeDerivs);

   double derivByParam(Arg &, const std::string &paramName)
      throw(ParameterNotFound);

   //! must re-implement here since this Statistic does not
   //! comprise individual Sources
   void setParams(std::vector<Parameter> &params) 
      throw(LikelihoodException, ParameterNotFound) {
      Function::setParams(params);
   }

private:

   double m_prefactor;

   void init();

};

} // namespace Likelihood

