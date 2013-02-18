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

#include <cmath>
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

class KahanAccumulator {
  // Kahan summation algorithm - summation with correction
  // Reference: http://en.wikipedia.org/wiki/Kahan_summation_algorithm

public:

  KahanAccumulator(double S=0, double C=0): m_S(S), m_C(C) { } 
  KahanAccumulator(const KahanAccumulator& o): m_S(o.m_S), m_C(o.m_C) { }

  KahanAccumulator & operator= (const KahanAccumulator& o) { 
     m_S=o.m_S; 
     m_C=o.m_C;
     return *this;
  }

  void reset(double s=0) { m_S = s; m_C = 0; }

  void add(double x)  {
     double Y(x-m_C);
     double T(m_S+Y);
     m_C=(T-m_S)-Y;
     m_S=T;
  }

  void add(const KahanAccumulator& x) {
    add(-x.m_C);
    add(x.m_S);
  }

  void addLog(double x) { add(std::log(x)); }

  void addLog(const KahanAccumulator& x) { 
     double log_sum = std::log(x.m_S);
     add(-x.m_C/x.m_S);
     add(log_sum);
  }
  
  double sum() const { return m_S; }
  double correction() const { return m_C; }

  // For compatibility with Accumulator
  double total() { double s = sum(); reset(); return s; }

  KahanAccumulator operator - () const { 
     KahanAccumulator a(-m_S, -m_C);
     return a;
  }

private:
  double m_S;
  double m_C;
};


} // namespace Likelihood

#endif // Likelihood_Accumulator_h
