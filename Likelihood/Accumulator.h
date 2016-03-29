/**
 * @file Accumulator.h
 * @brief Partition data by size in logrithmic bins and accumulate in
 * bins before summing.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/Accumulator.h,v 1.3 2012/02/07 00:24:27 jchiang Exp $
 */

#ifndef Likelihood_Accumulator_h
#define Likelihood_Accumulator_h

#include <vector>

namespace Likelihood {

class Kahan_Accumulator {

public:
  Kahan_Accumulator()
    :m_sum(0.),
     m_c(0.){
  }

  ~Kahan_Accumulator(){;}

  // Add a value to the running sum
  void add(double value) {
    // The next line pick up the rounding error 
    double y = value - m_c;
    // This adds the value and the rounding error to a temp variable
    double t = m_sum + y; 
    // This next line caluclates the rounding error
    // Note that the parentheses are there to force the evalutation order
    m_c = (t - m_sum) - y;
    // This line moves the temp to the running sum
    m_sum = t;
  }
  
  // Get the total and reset the running sum
  double total() {
    double t = m_sum;
    m_sum = 0.;
    m_c = 0.;
    return t;
  }

private:

  double m_sum; // This is the running sum
  double m_c;   // This tracks the rounding error

};




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

   double m_max_xrange;
   
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
