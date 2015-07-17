/** 
 * @file HistND.h
 * @brief N-dimensional histogram.
*/
#include <stdexcept>
#include <iostream>

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


void HistND::setSlice(unsigned int idim,
		      const std::vector<unsigned int>& ivalues,
		      const std::vector<float>& data) {       
  if ( idim >= m_ndims ) {
    throw std::length_error("HistND::setSlice:\n"
			    "idim >= histogram dimension.");
  }
  if ( ivalues.size() != m_ndims ) {
    throw std::length_error("HistND::setSlice:\n"
                            "length of ivalues != histogram dimension ");
  }
  unsigned int dimSize(0);
  if ( idim == m_ndims -1 ) {
    dimSize = m_data.size() / m_strides[idim];
  } else {
    dimSize = m_strides[idim+1] / m_strides[idim];
  }
  if ( data.size() != dimSize ) {
    throw std::length_error("HistND::setSlice:\n"
                            "length of input data != histogram slice size.");
  }
  // ivalues_start is a vector of all the indices for the first element in the slice
  std::vector<long> ivalues_start(ivalues.size(),0);
  std::copy(ivalues.begin(), ivalues.end(), ivalues_start.begin());
  ivalues_start[idim] = 0;

  // Get the starting bin and the strides, then just move along the strides and fill m_data
  long idx = binIndex(ivalues_start);
  long stride = m_strides[idim];
  for ( unsigned int i(0); i < data.size(); i++, idx+= stride ) {
    m_data[idx] = data[i];
  }

}

void HistND::setSlice(unsigned int idim,
		      const std::vector<unsigned int>& ivalues,
		      const std::vector<double>& data) {       
  if ( idim >= m_ndims ) {
    throw std::length_error("HistND::setSlice:\n"
			    "idim >= histogram dimension.");
  }
  if ( ivalues.size() != m_ndims ) {
    throw std::length_error("HistND::setSlice:\n"
                            "length of ivalues != histogram dimension ");
  }
  unsigned int dimSize(0);
  if ( idim == m_ndims -1 ) {
    dimSize = m_data.size() / m_strides[idim];
  } else {
    dimSize = m_strides[idim+1] / m_strides[idim];
  }
  if ( data.size() != dimSize ) {
    throw std::length_error("HistND::setSlice:\n"
                            "length of input data != histogram slice size.");
  }
  // ivalues_start is a vector of all the indices for the first element in the slice
  std::vector<long> ivalues_start(ivalues.size(),0);
  std::copy(ivalues.begin(), ivalues.end(), ivalues_start.begin());
  ivalues_start[idim] = 0;

  // Get the starting bin and the strides, then just move along the strides and fill m_data
  long idx = binIndex(ivalues_start);
  long stride = m_strides[idim];
  for ( unsigned int i(0); i < data.size(); i++, idx+= stride ) {
    m_data[idx] = data[i];
  }

}

void HistND::getSlice(unsigned int idim,
		      const std::vector<unsigned int>& ivalues,
		      std::vector<float>& data) {
  if ( idim >= m_ndims ) {
    throw std::length_error("HistND::getSlice:\n"
                            "idim >= histogram dimension.");
  }
  if ( ivalues.size() != m_ndims ) {
    throw std::length_error("HistND::getSlice:\n"
                            "length of ivalues != histogram dimension");
  }
  unsigned int dimSize(0);
  if ( idim == m_ndims -1 ) {
    dimSize = m_data.size() / m_strides[idim];
  } else {
    dimSize = m_strides[idim+1] / m_strides[idim];
  }
  data.clear();
  data.reserve(dimSize);

  // ivalues_start is a vector of all the indices for the first element in the slice 
  std::vector<long> ivalues_start(ivalues.size(),0);
  std::copy(ivalues.begin(), ivalues.end(), ivalues_start.begin());
  ivalues_start[idim] = 0;

  // Get the starting bin and the strides, then just move along the strides and fill the data vector
  long idx = binIndex(ivalues_start);
  long stride = m_strides[idim];
  for ( unsigned int i(0); i < dimSize; i++, idx+= stride ) {
    data.push_back(m_data[idx]);
  }

}

void HistND::getSlice(unsigned int idim,
		      const std::vector<unsigned int>& ivalues,
		      std::vector<double>& data) {
  if ( idim >= m_ndims ) {
    throw std::length_error("HistND::getSlice:\n"
                            "idim >= histogram dimension.");
  }
  if ( ivalues.size() != m_ndims ) {
    throw std::length_error("HistND::getSlice:\n"
                            "length of ivalues != histogram dimension");
  }
  unsigned int dimSize(0);
  if ( idim == m_ndims -1 ) {
    dimSize = m_data.size() / m_strides[idim];
  } else {
    dimSize = m_strides[idim+1] / m_strides[idim];
  }
  data.clear();
  data.reserve(dimSize);

  // ivalues_start is a vector of all the indices for the first element in the slice 
  std::vector<long> ivalues_start(ivalues.size(),0);
  std::copy(ivalues.begin(), ivalues.end(), ivalues_start.begin());
  ivalues_start[idim] = 0;

  // Get the starting bin and the strides, then just move along the strides and fill the data vector
  long idx = binIndex(ivalues_start);
  long stride = m_strides[idim];
  for ( unsigned int i(0); i < dimSize; i++, idx+= stride ) {
    data.push_back(m_data[idx]);
  }

}

HistND* HistND::sumRange(unsigned int idim, unsigned int firstBin, unsigned int lastBin) const {
  if ( idim >= m_ndims ) {
    throw std::length_error("HistND::sumRange:\n"
                            "idim >= histogram dimension.");
    return 0;
  }
  if ( firstBin >= lastBin ) {
    throw std::length_error("HistND::sumRange:\n"
                            "firstBin >= lastBin");
    return 0;
  }  
  unsigned int dimSize(0);
  if ( idim == m_ndims -1 ) {
    dimSize = m_data.size() / m_strides[idim];
  } else {
    dimSize = m_strides[idim+1] / m_strides[idim];
  }
  if ( firstBin >= dimSize ) {
    throw std::length_error("HistND::sumRange:\n"
                            "firstBin >= histogram axis size.");
    return 0;
  }  
  if ( lastBin > dimSize ) {
    throw std::length_error("HistND::sumRange:\n"
                            "lastBin > histogram axis size.");
    return 0;
  }  

  std::vector<evtbin::Binner*> binners;
  for ( unsigned int i(0); i < m_ndims; i++ ) {
    if ( i == idim ) continue;
    binners.push_back(const_cast<evtbin::Binner*>(getBinners()[i]));
  }
  
  HistND* retHist = new HistND(binners);
  std::vector<float> outData(retHist->data().size(),0.);
  
  unsigned int stride = m_strides[idim];
  for ( unsigned int iout(0); iout < outData.size(); iout++ ) {
    unsigned int idx = iout;
    for ( unsigned int j(firstBin); j < lastBin; j++, idx+= stride) {
      outData[iout] = m_data[idx];
    }
  }
  retHist->setData(outData);
  return retHist;
}


}
