/**
 * @file RunParams.h
 * @brief Wrapper class for the hoops interface.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/RunParams.h,v 1.7 2003/11/26 01:51:03 jchiang Exp $
 */

#ifndef Likelihood_RunParams_h
#define Likelihood_RunParams_h

#include "hoops/hoops.h"
#include "hoops/hoops_limits.h"
#include "hoops/hoops_pil_factory.h"
#include "hoops/hoops_pil.h"
#include "hoops/hoops_exception.h"

#include <map>
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/RunParams.h,v 1.7 2003/11/26 01:51:03 jchiang Exp $
 */

class RunParams {

public:

   RunParams() : m_prompter(0) {}

   RunParams(int iargc, char* argv[]);

   ~RunParams();

   template <typename T> 
   void getParam(const std::string &name, T &value) const {
// Failure to return an object of the requested type T is a logic
// error in the calling routine and is therefore considered
// unrecoverable in this context, hence the assertions.
      try {
         T my_value = (m_prompter->Group())[name];
         value = my_value;
      } catch (hoops::Hexception &eObj) {
         std::cout << "HOOPS exception: " << eObj.Msg() << "\n"
                   << "Code: " << eObj.Code() << std::endl;
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

   std::map<int, std::string> m_pil_errors;

   void errorCodes();

};

} // namespace Likeliood

#endif // Likelihood_RunParams_h
