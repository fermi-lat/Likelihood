/** 
 * @file Parameter.h
 * @brief Declaration of Parameter and OutOfBounds classes
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Parameter.h,v 1.11 2003/05/21 23:11:07 jchiang Exp $
 */
#pragma warning(disable:4290)
#ifndef Parameter_h
#define Parameter_h

#include <vector>
#include <string>
#include <cmath>
#include "Likelihood/LikelihoodException.h"

namespace Likelihood {

/** 
 * @class Parameter
 *
 * @brief Model parameters are identified by a name with flags to
 * indicate if it's free and with upper and lower bounds.
 *
 * ToDo: Need to check optimizer implementations for unbounded
 * parameters.  Does one simply set bounds at extremal (HUGE) values?
 * This may be a dangerous default behavior since some optimizers like
 * to evaluate at the bounds. Or can one simply set a flag to indicate
 * an unbounded parameter?
 *
 * @authors J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Parameter.h,v 1.11 2003/05/21 23:11:07 jchiang Exp $ 
 */

class OutOfBounds;

class Parameter {
    
public:
   
   Parameter() {init(std::string(""), 0., -HUGE, HUGE, true);};
   Parameter(const std::string &paramName, double paramValue, 
             bool isFree = true)
      {init(paramName, paramValue, -HUGE, HUGE, isFree);};
   Parameter(const std::string &paramName, double paramValue, 
             double minValue, double maxValue, bool isFree = true)
      {init(paramName, paramValue, minValue, maxValue, isFree);};

// need only member-wise copying
//   Parameter(const Parameter&);

   ~Parameter(){}

   //! name access
   void setName(const std::string &paramName) {m_name = paramName;};
   std::string getName() const {return m_name;}
   
   //! value access
   void setValue(double value) throw(OutOfBounds);
   double getValue() const {return m_value;}

   //! scale access
   void setScale(double scale) {m_scale = scale;}
   double getScale() const {return m_scale;}

   //! "true" value access
   void setTrueValue(double trueValue) throw(OutOfBounds);
   double getTrueValue() const {return m_value*m_scale;}

   //! bounds access
   void setBounds(double minValue, double maxValue) throw(OutOfBounds);
   void setBounds(const std::pair<double, double> &boundValues) 
      throw(OutOfBounds) {setBounds(boundValues.first, boundValues.second);}
   std::pair<double, double> getBounds();

   //! free flag access
   void setFree(bool free) {m_free = free;};
   bool isFree() const {return m_free;};

private:

   //! set all the Parameter values (with default scaling)
   void init(const std::string &paramName, double paramValue, 
             double minValue, double maxValue, bool isFree = true)
      {m_name = paramName; m_value = paramValue; m_minValue = minValue;
      m_maxValue = maxValue; m_free = isFree; m_scale = 1;}

   //! parameter name
   std::string m_name;

   //! its value
   double m_value;

   //! its scale factor
   double m_scale;

   //! lower bound
   double m_minValue;

   //! upper bound
   double m_maxValue;

   //! flag to indicate free or fixed
   bool m_free;

};

/**
 * @class OutOfBounds
 *
 * @brief Exception class to ensure set[True]Value and setBounds methods
 * behave consistently with regard to existing values.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Parameter.h,v 1.11 2003/05/21 23:11:07 jchiang Exp $
 */

class OutOfBounds : public LikelihoodException {

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

#endif // Parameter_h
