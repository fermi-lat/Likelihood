/**
 * @file AppHelpers.h
 * @brief Class of "helper" methods for the Likelihood applications.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/AppHelpers.h,v 1.8 2004/12/08 00:31:11 jchiang Exp $
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/AppHelpers.h,v 1.8 2004/12/08 00:31:11 jchiang Exp $
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
   void setRoi(const std::string & filename="", const std::string & ext="",
               bool strict=true);

   static void checkOutputFile(bool clobber, const std::string & filename);

   void checkOutputFile() {
      checkOutputFile(m_pars["clobber"], m_pars["outfile"]);
   }

   const std::vector<std::string> & scFiles() const {return m_scFiles;}

   static void checkCuts(const std::string & file1, const std::string ext1,
                         const std::string & file2, const std::string ext2);

protected:

   st_app::AppParGroup & m_pars;
   optimizers::FunctionFactory * m_funcFactory;
   std::vector<std::string> m_scFiles;

   void prepareFunctionFactory();
   void createResponseFuncs();
};

} // namespace Likelihood

#endif // Likelihood_AppHelpers
