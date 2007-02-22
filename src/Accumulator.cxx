/**
 * @file Accumulator.cxx
 * @brief Partition data by size in logrithmic bins and accumulate in
 * bins before summing.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Accumulator.cxx,v 1.1 2007/02/19 18:05:32 jchiang Exp $
 */

#include <cmath>

#include <algorithm>

#include "Likelihood/Accumulator.h"

namespace Likelihood {

Accumulator::Accumulator() : m_setOuterBounds(true), m_first(true), 
                             m_total(0), m_npts(10) {}

void Accumulator::add(double value) {
   double absvalue(std::abs(value));
   if (m_setOuterBounds && absvalue !=0) {
      m_xmin = absvalue;
      m_xmax = absvalue;
      m_setOuterBounds = false;
      m_total = 0;
   } else {
      if (absvalue < m_xmin && absvalue !=0) {
         m_xmin = absvalue;
      }
      if (absvalue > m_xmax) {
         m_xmax = absvalue;
      }
   }
   if (m_first) {
      m_total += value;
   } else {
      m_partialSums.at(indx(absvalue)).add(value);
   }
}

double Accumulator::total() {
   if (m_first) {
      m_first = false;
      m_partialSums.resize(m_npts);
      return m_total;
   }
   std::vector<double> sums;
   for (size_t i(0); i < m_partialSums.size(); i++) {
      sums.push_back(m_partialSums.at(i).total());
   }
   std::stable_sort(sums.begin(), sums.end());
   m_total = 0;
   for (size_t i(0); i < sums.size(); i++) {
      m_total += sums.at(i);
   }
   m_partialSums.clear();
   m_partialSums.resize(m_npts);
   return m_total;
}

size_t Accumulator::indx(double xx) {
   if (xx <= m_xmin) {
      return 0;
   } else if (xx >= m_xmax) {
      return m_partialSums.size() - 1;
   }
   double xstep(std::log(m_xmax/m_xmin)/(m_npts - 1));
   return static_cast<size_t>(std::log(xx/m_xmin)/xstep);
}

} // namespace Likelihood
