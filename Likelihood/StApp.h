/**
 * @file StApp.h
 * @brief Template for st_app::IApp.
 * $Header$
 */

#ifndef Likelihood_StApp_h
#define Likelihood_StApp_h

#include "st_app/IApp.h"
#include "hoops/hoops_prompt_group.h"

namespace Likelihood {

/**
 * @class StApp
 * @brief Class template of boiler-plate code expected by st_app. 
 * 
 * $Header$
 */
template<typename T>
class StApp : public st_app::IApp {
public:
   StApp(const std::string &appName) : st_app::IApp(appName) {}
   virtual void run() {
      hoopsPrompt();
      hoopsSave();
      hoops::IParGroup & pars = hoopsGetParGroup();
      T my_object(pars);
      my_object.run();
   }
};

} // namespace Likelihood

#endif // Likelihood_StApp_h
