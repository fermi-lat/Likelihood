/**
 * @file AppBase.h
 * @brief Base class for Likelihood applications providing common
 * functionality.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef Likelihood_AppBase
#define Likelihood_AppBase

#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <string>

#include "st_app/IApp.h"

#include "optimizers/FunctionFactory.h"

namespace Likelihood {

/**
 * @class AppBase
 * @brief A subclass of st_app::IApp that serves as a base
 * class for Likelihood applications.
 *
 * @author J. Chiang
 *
 * $Header$
 */

class AppBase : public st_app::IApp {

public:

   virtual ~AppBase() throw() {}

   virtual void setUp();
   virtual void tearDown();

protected:

   optimizers::FunctionFactory m_funcFactory;

   AppBase(const std::string & appName) : st_app::IApp(appName) {}

   virtual void prepareFunctionFactory();
   virtual void setRoi();
   virtual void readScData();
   virtual void createResponseFuncs();
   virtual void readExposureMap();

};

} // namespace Likelihood

#endif // Likelihood_AppBase
