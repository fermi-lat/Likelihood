/**
 * @file AppHelpers.h
 * @brief Class of "helper" methods for the Likelihood applications.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/AppHelpers.h,v 1.1 2004/04/21 20:58:38 jchiang Exp $
 */

#ifndef Likelihood_AppHelpers
#define Likelihood_AppHelpers

#include <string>

#include "st_app/AppParGroup.h"

#include "optimizers/FunctionFactory.h"

namespace Likelihood {

/**
 * @class AppHelpers
 * @brief The methods in this class calls various static methods for
 * reading in spacecraft data, defining the region-of-interest, and
 * preparing the response functions --- all standard tasks which must
 * be performed as part of any Likelihood analysis.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/AppHelpers.h,v 1.1 2004/04/21 20:58:38 jchiang Exp $
 */

class AppHelpers {

public:

   AppHelpers(st_app::AppParGroup & pars) : m_pars(pars),
      m_funcFactory(0) {
      prepareFunctionFactory();
      readScData();
      createResponseFuncs();
   }

   ~AppHelpers() {
      delete m_funcFactory;
   }

   optimizers::FunctionFactory & funcFactory();

   void readExposureMap();
   void setRoi();

protected:

   st_app::AppParGroup & m_pars;
   optimizers::FunctionFactory * m_funcFactory;

   void prepareFunctionFactory();
   void readScData();
   void createResponseFuncs();
};

} // namespace Likelihood

#endif // Likelihood_AppHelpers
