/**
 * @file RunParams.h
 * @brief Wrapper class for the hoops interface.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef Likelihood_RunParams_h
#define Likelihood_RunParams_h

#include "hoops/hoops.h"
#include "hoops/hoops_limits.h"
#include "hoops/hoops_pil_factory.h"
#include "hoops/hoops_pil.h"
#include "hoops/hoops_exception.h"

#include <string>

namespace Likelihood {

/**
 * @class RunParams
 * @brief This class wraps the HOOPS interface and retrieves command-line
 * parameters for the standalone likelihood application.
 *
 * @author J. Chiang
 *
 * $Header$
 */

class RunParams {

public:

   RunParams(int iargc, char* argv[]);

   ~RunParams();

//    template <typename T> T operator()(const std::string &name) const {
//       hoops::IParGroup & pg = m_prompter->Group();
//       T value = pg[name];
//       return value;
//    }
   
   std::string string_par(const std::string &name) const {
      try {
         std::string my_string = (m_prompter->Group())[name];
         return my_string;
      } catch (...) {
         std::cout << name << std::endl;
         assert(false);
      }
   }

   double double_par(const std::string &name) const {
      try {
         double value = (m_prompter->Group())[name];
         return value;
      } catch (...) {
         std::cout << name << std::endl;
         assert(false);
      }
   }

   int int_par(const std::string &name) const {
      try {
         int value = (m_prompter->Group())[name];
         return value;
      } catch (hoops::Hexception &eObj) {
         std::cout << eObj.Msg() << "\n"
                   << eObj.Code() << std::endl;
         assert(false);
      } catch (...) {
         std::cout << name << std::endl;
         assert(false);
      }
   }

   long long_par(const std::string &name) const {
      try {
         long value = (m_prompter->Group())[name];
         return value;
      } catch (hoops::Hexception &eObj) {
         std::cout << eObj.Msg() << "\n"
                   << eObj.Code() << std::endl;
         assert(false);
      } catch (...) {
         std::cout << name << std::endl;
         assert(false);
      }
   }

private:

   hoops::IParPrompt * m_prompter;

};

} // namespace Likeliood

#endif // Likelihood_RunParams_h
