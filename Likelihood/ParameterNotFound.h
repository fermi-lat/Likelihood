/**
 * @file ParameterNotFound.h
 * @brief Declaration and definition of ParameterNotFound exception class.
 * @author J. Chiang
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ParameterNotFound.h,v 1.1 2003/06/19 20:16:04 jchiang Exp $
 */

#ifndef Likelihood_ParameterNotFound_h
#define Likelihood_ParameterNotFound_h

#include <sstream>
#include "Likelihood/LikelihoodException.h"

namespace Likelihood {

/**
 * @class ParameterNotFound
 *
 * @brief A class that returns a standard error message for
 * Parameters looked for but not found in the desired Function.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ParameterNotFound.h,v 1.1 2003/06/19 20:16:04 jchiang Exp $
 */

class ParameterNotFound : public LikelihoodException {

public:
   ParameterNotFound(const std::string &paramName, 
                     const std::string &funcName,
                     const std::string &routineName) {
      std::ostringstream errorMessage;
      errorMessage << "Function::" << routineName << ": \n"
                   << "A Parameter named " << paramName
                   << " is not a Parameter of Function "
                   << funcName << "\n";
      m_what = errorMessage.str();
      m_code = 0;
   }
};

} // namespace Likelihood

#endif // Likelihood_ParameterNotFound_h
