/** @file Parameter.cxx
 * @brief Parameter class implementation
 *
 * $Header:
 */

#include <vector>
#include <string>

#include "../Likelihood/Parameter.h"

namespace Likelihood {

//! copy constructor
Parameter::Parameter(const Parameter &param) {
   m_name = param.m_name;
   m_value = param.m_value;
   m_minValue = param.m_minValue;
   m_maxValue = param.m_maxValue;
   m_free = param.m_free;
}

//! return bounds as a pair
std::pair<double, double> Parameter::getBounds() {
   std::pair<double, double> my_Bounds(m_minValue, m_maxValue);
   return my_Bounds;
}

} // namespace Likelihood
