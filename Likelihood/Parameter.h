/** @file Parameter.h
 * @brief Declaration of Parameter class
 * $Header:$
 */

#ifndef Parameter_h
#define Parameter_h

#include <vector>
#include <string>
#include <cmath>

namespace Likelihood {

/** 
 * @class Parameter
 *
 * @brief Model parameters are identified by a name with flags to
 * indicate if its free and with upper and lower bounds.
 *
 * ToDo: Need to check optimizer implementations for unbounded
 * parameters.  Does one simply set bounds at extremal (HUGE)
 * values?  This may be a dangerous default behavior since some
 * optimizers like to evaluate at the bounds. Or can one simply set a
 * flag to indicate an unbounded parameter?
 *
 * @authors J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Parameter.h,v 1.4 2003/03/11 22:01:15 jchiang Exp $ */

class Parameter {
    
public:
   
   Parameter() {m_init(std::string(""), 0., -HUGE, HUGE, true);};
   Parameter(const std::string &paramName, double paramValue, 
             bool isFree = true)
      {m_init(paramName, paramValue, -HUGE, HUGE, isFree);};
   Parameter(const std::string &paramName, double paramValue, 
             double minValue, double maxValue, bool isFree = true)
      {m_init(paramName, paramValue, minValue, maxValue, isFree);};
   Parameter(const Parameter&);
   ~Parameter(){}

   //! name access
   void setName(const std::string &paramName) {m_name = paramName;};
   std::string getName() const {return m_name;}
   
   //! value access
   void setValue(double value) {m_value = value;}
   double getValue() const {return m_value;}

   //! bounds access
   void setBounds(double minValue, double maxValue)
      {m_minValue = minValue; m_maxValue = maxValue;};
   void setBounds(const std::pair<double, double> &boundValues)
      {setBounds(boundValues.first, boundValues.second);};
   std::pair<double, double> getBounds();

   //! free flag access
   void setFree(bool free) {m_free = free;};
   bool isFree() const {return m_free;};

private:

   //! set all the Parameter values
   void m_init(const std::string &paramName, double paramValue, 
               double minValue, double maxValue, bool isFree = true)
      {m_name = paramName; m_value = paramValue; m_minValue = minValue;
      m_maxValue = maxValue; m_free = isFree;}

   //! parameter name
   std::string m_name;

   //! its value
   double m_value;

   //! lower bound
   double m_minValue;

   //! upper bound
   double m_maxValue;

   //! flag to indicate free or fixed
   bool m_free;

};

} // namespace Likelihood

#endif // Parameter_h
