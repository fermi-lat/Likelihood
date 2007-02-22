/**
 * @file Accumulator.h
 * @brief Partition data by size in logrithmic bins and accumulate in
 * bins before summing.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef Likelihood_Accumulator_h
#define Likelihood_Accumulator_h

#include <vector>

namespace Likelihood {

class Accumulator {

public:

   Accumulator();

   void add(double value);

   double total();

private:

   bool m_setOuterBounds;
   bool m_first;
   
   double m_xmin;
   double m_xmax;
   double m_total;

   size_t m_npts;
   
   /**
    * @class SumByNumber
    * @brief Sum in standard chunk sizes.
    */
   class SumByNumber {
   public:
      SumByNumber(size_t num=100) : m_num(num), m_counter(0), m_chunks(1, 0) {}
      void add(double value) {
         m_chunks.back() += value;
         m_counter++;
         if (m_counter >= m_num) {
            m_chunks.push_back(0);
            m_counter = 0;
         }
      }
      double total() {
         double my_total(0);
         for (size_t i(0); i < m_chunks.size(); i++) {
            my_total += m_chunks.at(i);
         }
         m_chunks.resize(1, 0);
         m_counter = 0;
         return my_total;
      }
   private:
      size_t m_num;
      size_t m_counter;
      std::vector<double> m_chunks;
   };

   std::vector<SumByNumber> m_partialSums;

   size_t indx(double xx);

};

} // namespace Likelihood

#endif // Likelihood_Accumulator_h
