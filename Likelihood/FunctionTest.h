/**
 * @file FunctionTest.h
 * @brief Declaration of unit test code for Function class hierarchy
 * @author J. Chiang
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/FunctionTest.h,v 1.2 2003/05/29 20:10:30 jchiang Exp $
 */


#ifndef FunctionTest_h
#define FunctionTest_h

#include "Likelihood/LikelihoodException.h"
#include "Likelihood/Parameter.h"
#include "Likelihood/Arg.h"
#include "Likelihood/Function.h"

namespace Likelihood {

/**
 * @class FunctionTest
 * @brief Unit test code for Function class hierarchy.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/FunctionTest.h,v 1.2 2003/05/29 20:10:30 jchiang Exp $
 */

class FunctionTest {

public:

   FunctionTest(Function &func, const std::string &name) {
      m_func = &func;
      m_func->setName(name);
      assert(m_func->getName() == name);
      m_func->getParams(m_originalParameters);
   }

   ~FunctionTest() {}

   void parameters(const std::vector<Parameter> &params) 
      throw(LikelihoodException);

   void freeParameters(const std::vector<Parameter> &params) 
      throw(LikelihoodException);

   void funcEvaluations(const std::vector<Arg*> &arguments,
                        const std::vector<double> &returnValues)
      throw(LikelihoodException);

   void derivatives(const std::vector<Arg*> &arguments,
                    double eps = 1e-5) throw(LikelihoodException);
      
private:

   Function *m_func;
   std::vector<Parameter> m_originalParameters;

};

} // namespace FunctionTest

#endif // FunctionTest_h
