/**
 * @file AppHelpers.h
 * @brief Class of "helper" methods for the Likelihood applications.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/AppHelpers.h,v 1.2 2004/06/05 15:22:14 jchiang Exp $
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/AppHelpers.h,v 1.2 2004/06/05 15:22:14 jchiang Exp $
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

   const std::vector<std::string> & scFiles() const {return m_scFiles;}

protected:

   st_app::AppParGroup & m_pars;
   optimizers::FunctionFactory * m_funcFactory;
   std::vector<std::string> m_scFiles;

   void prepareFunctionFactory();
   void readScData();
   void createResponseFuncs();
};

} // namespace Likelihood

#endif // Likelihood_AppHelpers
