/** \file HistND.h
    \brief Three dimensional histogram.
*/
#ifndef Likelihood_HistND_h
#define Likelihood_HistND_h

#include <vector>

#include "evtbin/Hist.h"

namespace evtbin {
   class Binner;
}

namespace Likelihood {

   /** \class HistND
       \brief N dimensional histogram.
   */
   class HistND : public evtbin::Hist {
   public:

      /** \brief Create a three dimensional histogram which uses the given 
                 binner objects to determine the indices.
          \param binners vector of references to binners.  The size of this 
                 vector sets the dimensionality of the histogram.
      */
      HistND(const std::vector<Binner *> & binners);

      HistND(const HistND & rhs);
      
      virtual ~HistND() throw();
     
      /** \brief Increment the bin appropriate for the given value.
                 This is generic for N-dimensional histograms.
          \param value Vector giving the value being binned. The vector 
                 must have at least as many values as the dimensionality 
                 of the histogram.
      */
      virtual void fillBin(const std::vector<double> & values, 
                           double weight = 1.);

      long binIndex(const std::vector<double> & values) const;

      long binIndex(const std::vector<long> & ivalues) const;

      double operator[](long index) const;

      const std::vector<double> & data() const {return m_data;}

      void setData(const std::vector<double> & data) {m_data = data;}

      HistND * clone() const {return new HistND(*this);}

   private:

      std::vector<double> m_data;
     
      unsigned int m_ndims;
     
      std::vector<unsigned int> m_strides;

   };

}

#endif
