/**
 * @file OptEM.h
 * @brief Expectation maximization for fitting Sources.
 * @author P. L. Nolan
 *
 * $Header$
 */

#ifndef Likelihood_OptEM_h
#define Likelihood_OptEM_h

#include "LogLike.h"

namespace Likelihood {

/**
 * @class OptEM
 * @brief Expectation maximization class.
 *
 * Although this class inherits from LogLike, note that it has its own
 * findMin(...) method and so should not be used as an argument to an
 * optimizer::Optimizer constructor.
 *
 * @author P. L. Nolan
 *
 * $Header$
 */

  class OptEM: public LogLike {
  public:
    OptEM() {};
    virtual ~OptEM() {};
    void findMin(const int verbose = 0);

  private:

  }; //class OptEM
} // namespace Likelihood

#endif // Likelihood_OptEM_h
