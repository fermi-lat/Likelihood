/**
 * @file OutOfBounds.h
 * @brief Declaration/definition of Parameter::OutOfBounds exception class
 * @author J. Chiang
 * $Header$
 */

#ifndef OutOfBounds_h
#define OutOfBounds_h

#include "LikelihoodException.h"

namespace Likelihood {

/**
 * @class OutOfBounds
 *
 * @brief Nested exception class to ensure set[True]Value and setBounds 
 * methods behave consistently with regard to existing values.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Parameter.h,v 1.18 2003/06/10 19:31:09 jchiang Exp $
 */

class Parameter::OutOfBounds : public LikelihoodException {

public:
   OutOfBounds(const std::string &errorString, double value, 
               double minValue, double maxValue, int code) : 
      LikelihoodException(errorString, code), m_value(value), 
      m_minValue(minValue), m_maxValue(maxValue) {}
   ~OutOfBounds() {}
   
   double value() {return m_value;}
   double minValue() {return m_minValue;}
   double maxValue() {return m_maxValue;}
   
   enum ERROR_CODES {VALUE_ERROR, BOUNDS_ERROR};

private:

   double m_value;
   double m_minValue;
   double m_maxValue;

};

} // namespace Likelihood

#endif // OutOfBounds_h
