/**
 * @file SparseVector.cxx
 * @brief Implmentation of a sparse vector (similar to an std::map, but with less functionality, and using less memory)
 * @author E. Charles
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/FileUtils.cxx,v 1.2 2016/09/21 22:42:40 echarles Exp $
 */


#include "Likelihood/SparseVector.h"

namespace Likelihood {

  static SparseVector<int> sparse_int;
  static SparseVector<float> sparse_float;
  static SparseVector<double> sparse_double;
 
} // namespace Likelihood
