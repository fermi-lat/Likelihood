/**
 * @file FunctionTest.h
 * @brief Declaration of unit test code for Function class hierarchy
 * @author J. Chiang
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/FunctionTest.h,v 1.6 2003/07/19 04:38:01 jchiang Exp $
 */


#ifndef Likelihood_FunctionTest_h
#define Likelihood_FunctionTest_h

#include "Likelihood/Parameter.h"
#include "Likelihood/Arg.h"
#include "Likelihood/Function.h"
#include "Likelihood/Exception.h"

namespace Likelihood {

/**
 * @class FunctionTest
 * @brief Unit test code for Function class hierarchy.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/FunctionTest.h,v 1.6 2003/07/19 04:38:01 jchiang Exp $
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
      throw(Exception);

   void freeParameters(const std::vector<Parameter> &params) 
      throw(Exception);

   void funcEvaluations(const std::vector<Arg*> &arguments,
                        const std::vector<double> &returnValues)
      throw(Exception);

   void derivatives(const std::vector<Arg*> &arguments,
                    double eps = 1e-5) throw(Exception);
      
private:

   Function *m_func;
   std::vector<Parameter> m_originalParameters;

};

} // namespace FunctionTest

#endif // Likelihood_FunctionTest_h
