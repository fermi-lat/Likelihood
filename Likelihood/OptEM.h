/**
 * @file OptEM.h
 * @brief Expectation maximization for fitting Sources.
 * @author P. L. Nolan
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/OptEM.h,v 1.5 2005/02/27 06:42:24 jchiang Exp $
 */

#ifndef Likelihood_OptEM_h
#define Likelihood_OptEM_h

#include "Likelihood/LogLike.h"
#include "Likelihood/Observation.h"

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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/OptEM.h,v 1.5 2005/02/27 06:42:24 jchiang Exp $
 */

  class OptEM: public LogLike {
  public:
    OptEM(const Observation & observation) : LogLike(observation) {}
    virtual ~OptEM() {}
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
