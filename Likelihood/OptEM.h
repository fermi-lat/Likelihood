/**
 * @file OptEM.h
 * @brief Expectation maximization for fitting Sources.
 * @author P. L. Nolan
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/OptEM.h,v 1.3 2003/11/22 00:06:12 pln Exp $
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/OptEM.h,v 1.3 2003/11/22 00:06:12 pln Exp $
 */

  class OptEM: public LogLike {
  public:
    OptEM() {};
    virtual ~OptEM() {};
    void findMin(const int verbose = 0);

  protected:
     virtual OptEM * clone() {
        return new OptEM(*this);
     }

  private:

  }; //class OptEM

  double chifunc(int);

} // namespace Likelihood

#endif // Likelihood_OptEM_h
