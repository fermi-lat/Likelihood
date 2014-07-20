/** 
 * @file HistND.h
 * @brief N-dimensional histogram.
*/
#include <stdexcept>

#include "evtbin/Binner.h"
#include "Likelihood/HistND.h"

namespace Likelihood {

HistND::HistND(const std::vector<evtbin::Binner *> & binners) {
// Copy binners and set size of the data array.
   m_ndims = binners.size();
   m_binners.resize(m_ndims);
   int num_bins(1);
   for (unsigned int i = 0; i < m_ndims; i++) {
      m_binners[i] = binners[i]->clone();
      num_bins *= binners[i]->getNumBins();
   }
   m_data.resize(num_bins, 0);

// Compute strides.
   m_strides.resize(m_ndims);
   m_strides[0] = 1;
   for (unsigned int i = 1; i < m_ndims; i++) {
      m_strides[i] = m_strides[i-1]*binners[i-1]->getNumBins();
   }
}

HistND::HistND(const HistND & rhs) : evtbin::Hist(rhs) {
   const evtbin::Hist::BinnerCont_t & rhs_binners = rhs.getBinners();
   long n_binners = rhs_binners.size();
   m_binners.resize(n_binners);
   for (long i = 0; i < n_binners; i++) {
      m_binners[i] = rhs_binners[i]->clone();
   }
   m_data = rhs.m_data;
   m_ndims = rhs.m_ndims;
   m_strides = rhs.m_strides;
}

HistND::~HistND() throw() {}

void HistND::fillBin(const std::vector<double> & values, double weight) {
   long indx = binIndex(values);
   if (indx >= 0) {
      m_data[indx] += weight;
   }
}

long HistND::binIndex(const std::vector<double> & values,
                      long border_size) const {
   if (values.size() != m_ndims) {
      throw std::length_error("HistND::binIndex:\n"
                              "Size of values vector does match "
                              "histogram dimension.");
   }
   std::vector<long> ivalues;
   for (unsigned int i = 0; i < values.size(); i++) {
      long bin_index = m_binners[i]->computeIndex(values[i]);
      if (bin_index < border_size || 
          bin_index >= m_binners[i]->getNumBins() - border_size) {
         return -1;
      } else {
         ivalues.push_back(bin_index);
      }
   }
   return binIndex(ivalues);
}

long HistND::binIndex(const std::vector<long> & ivalues) const {
   if (ivalues.size() != m_ndims) {
      throw std::length_error("HistND::binIndex:\n"
                              "Size of ivalues vector does match "
                              "histogram dimension.");
   }
   long indx(0);
   for (unsigned int i = 0; i < ivalues.size(); i++) {
      indx += m_strides[i]*ivalues[i];
   }
   return indx;
}

double HistND::operator[](long index) const {
   return m_data[index];
}

}
