/** 
 * @file Parameter.cxx
 * @brief Parameter class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Parameter.cxx,v 1.9 2003/06/10 23:58:51 jchiang Exp $
 */

#include <vector>
#include <string>

#include "Likelihood/Parameter.h"

namespace Likelihood {

void Parameter::setValue(double value) throw(OutOfBounds) {
   if (value >= m_minValue && value <= m_maxValue) {
      m_value = value;
   } else {
      throw OutOfBounds(
         "Attempt to set the value outside of existing bounds.", 
         value, m_minValue, m_maxValue, 
         static_cast<int>(OutOfBounds::VALUE_ERROR));
   }
}

void Parameter::setTrueValue(double trueValue) throw(OutOfBounds) {
   double value = trueValue/m_scale;
   setValue(value);
}

void Parameter::setBounds(double minValue, double maxValue) 
   throw(OutOfBounds) {
   if (m_value >= minValue && m_value <= maxValue) {
      m_minValue = minValue;
      m_maxValue = maxValue;
   } else {
      throw OutOfBounds(
         "Attempt to set bounds that exclude the existing value.", 
         m_value, minValue, maxValue, 
         static_cast<int>(OutOfBounds::BOUNDS_ERROR));
   }
}

// return bounds as a pair
std::pair<double, double> Parameter::getBounds() const {
   std::pair<double, double> my_Bounds(m_minValue, m_maxValue);
   return my_Bounds;
}

} // namespace Likelihood
