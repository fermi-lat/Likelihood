/** 
 * @file Parameter.h
 * @brief Declaration of Parameter and OutOfBounds classes
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Parameter.h,v 1.16 2003/06/10 18:18:29 burnett Exp $
 */

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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Parameter.h,v 1.16 2003/06/10 18:18:29 burnett Exp $ 
 */

class Parameter {

   class OutOfBounds;
    
public:
   
   Parameter() {init(std::string(""), 0., -HUGE, HUGE, true);}

   /**
    * @param paramName The name of the Parameter
    * @param paramValue The (scaled) value of the Parameter
    * @param minValue Parameter value lower bound
    * @param maxValue Parameter value upper bound
    * @param isFree true if the Parameter value is allowed to vary in a fit
    */
   Parameter(const std::string &paramName, double paramValue, 
             double minValue, double maxValue, bool isFree = true)
      {init(paramName, paramValue, minValue, maxValue, isFree);}

   Parameter(const std::string &paramName, double paramValue, 
             bool isFree = true)
      {init(paramName, paramValue, -HUGE, HUGE, isFree);}

// need only member-wise copying
//   Parameter(const Parameter&);

   ~Parameter(){}

   //! name access
   void setName(const std::string &paramName) {m_name = paramName;}
   std::string getName() const {return m_name;}
   
   //! value access
   void setValue(double value) throw(OutOfBounds);
   double getValue() const {return m_value;}

   /**
    * @brief Scale access.
    *
    * The true value of the Parameter is used in the Function
    * calculation.  Only the (apparent) value is intended to
    * accessible through the value accessor methods of the Function
    * class.
    */
   void setScale(double scale) {m_scale = scale;}
   double getScale() const {return m_scale;}

   //! "true" value access
   void setTrueValue(double trueValue) throw(OutOfBounds);
   double getTrueValue() const {return m_value*m_scale;}

   //! bounds access
   void setBounds(double minValue, double maxValue) throw(OutOfBounds);
   void setBounds(const std::pair<double, double> &boundValues) 
      throw(OutOfBounds) {setBounds(boundValues.first, boundValues.second);}
   std::pair<double, double> getBounds() const;

   //! free flag access
   void setFree(bool free) {m_free = free;}
   bool isFree() const {return m_free;}

/**
 * @class OutOfBounds
 *
 * @brief Nested exception class to ensure set[True]Value and setBounds 
 * methods behave consistently with regard to existing values.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Parameter.h,v 1.16 2003/06/10 18:18:29 burnett Exp $
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


} // namespace Likelihood

#endif // Parameter_h
