/**
 * @file ParameterNotFound.h
 * @brief Declaration and definition of Function::ParameterNotFound 
 * exception class.
 * @author J. Chiang
 * $Header$
 */

#ifndef ParameterNotFound_h
#define ParameterNotFound_h

#include "LikelihoodException.h"

namespace Likelihood {

/**
 * @class Function::ParameterNotFound
 *
 * @brief Nested class that returns a standard error message for
 * Parameters looked for but not found in the desired Function.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Function.h,v 1.26 2003/06/10 20:20:42 jchiang Exp $
 */

class Function::ParameterNotFound : public LikelihoodException {

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

#endif // ParameterNotFound_h
