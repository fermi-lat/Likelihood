/**
 * @file RunParams.h
 * @brief Wrapper class for the hoops interface.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/RunParams.h,v 1.4 2003/11/10 12:54:39 cohen Exp $
 */

#ifndef Likelihood_RunParams_h
#define Likelihood_RunParams_h

#include "hoops/hoops.h"
#include "hoops/hoops_limits.h"
#include "hoops/hoops_pil_factory.h"
#include "hoops/hoops_pil.h"
#include "hoops/hoops_exception.h"

#include <string>
#include <iostream>

namespace Likelihood {

/**
 * @class RunParams
 * @brief This class wraps the HOOPS interface and retrieves command-line
 * parameters for the standalone likelihood application.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/RunParams.h,v 1.4 2003/11/10 12:54:39 cohen Exp $
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
      } catch (hoops::Hexception &eObj) {
         std::cout << "HOOPS exception: " << eObj.Msg() << "\n"
                   << "Code: " << eObj.Code() << std::endl;
         assert(false);
      } catch (...) {
         std::cout << name << std::endl;
         assert(false);
      }
   }

   double double_par(const std::string &name) const {
      try {
         double value = (m_prompter->Group())[name];
         return value;
      } catch (hoops::Hexception &eObj) {
         std::cout << "HOOPS exception: " << eObj.Msg() << "\n"
                   << "Code: " << eObj.Code() << std::endl;
         assert(false);
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
         std::cout << "HOOPS exception: " << eObj.Msg() << "\n"
                   << "Code: " << eObj.Code() << std::endl;
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

   bool bool_par(const std::string &name) const {
      try {
         bool value = (m_prompter->Group())[name];
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

// Static methods for parsing filename parameters.

   /// @param filename This can be a FITS file or a list of FITS
   ///                 files.  
   /// @param files On return, this contains the list of FITS files,
   ///              or the filename itself as its only element.
   static void resolve_fits_files(std::string filename, 
                                  std::vector<std::string> &files);

   /// @param inputFile The name of the file containing a list of strings
   ///                  separated by '\n'.
   /// @param lines On return, this contains the list of strings.
   static void readLines(std::string inputFile, 
                         std::vector<std::string> &lines);

private:

   hoops::IParPrompt * m_prompter;

};

} // namespace Likeliood

#endif // Likelihood_RunParams_h
