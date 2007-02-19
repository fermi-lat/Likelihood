/**
 * @file Accumulator.h
 * @brief Partition data by size in logrithmic bins and accumulate in
 * bins before summing.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef Accumulator_h
#define Accumulator_h

#include <vector>

namespace Likelihood {

class Accumulator {

public:

   Accumulator(size_t npts=10) : m_first(true), m_setOuterBounds(true),
      m_npts(npts) {}

   void add(double value);

   double total();
   
private:
   
   bool m_first;
   bool m_setOuterBounds;
   size_t m_npts;
   double m_xmin;
   double m_xmax;
   double m_total;
   std::vector<double> m_partialSums;

   size_t indx(double xx);
};

} // namespace Likelihood

#endif // Accumulator_h
