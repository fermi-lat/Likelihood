/**
 * @file AppHelpers.h
 * @brief Class of "helper" methods for the Likelihood applications.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/AppHelpers.h,v 1.6 2004/11/01 06:27:38 jchiang Exp $
 */

#ifndef Likelihood_AppHelpers
#define Likelihood_AppHelpers

#include <string>

#include "st_app/AppParGroup.h"

#include "optimizers/FunctionFactory.h"

namespace Likelihood {

/**
 * @class AppHelpers
 * @brief The methods in this class call various static methods for
 * reading in spacecraft data, defining the region-of-interest, and
 * preparing the response functions --- all standard tasks which must
 * be performed as part of any Likelihood analysis.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/AppHelpers.h,v 1.6 2004/11/01 06:27:38 jchiang Exp $
 */

class AppHelpers {

public:

   AppHelpers(st_app::AppParGroup & pars);

   ~AppHelpers() {
      delete m_funcFactory;
   }

   optimizers::FunctionFactory & funcFactory();

   void readScData();
   void readExposureMap();
   void setRoi();

   static void checkOutputFile(bool clobber, const std::string & filename);

   void checkOutputFile() {
      checkOutputFile(m_pars["clobber"], m_pars["outfile"]);
   }

   const std::vector<std::string> & scFiles() const {return m_scFiles;}

protected:

   st_app::AppParGroup & m_pars;
   optimizers::FunctionFactory * m_funcFactory;
   std::vector<std::string> m_scFiles;

   void prepareFunctionFactory();
   void createResponseFuncs();
};

} // namespace Likelihood

#endif // Likelihood_AppHelpers
