#include "LogLike.h"

namespace Likelihood {

  class OptEM: public LogLike {
  public:
    OptEM() {};
    virtual ~OptEM() {};
    void findMin(const int verbose = 0);

  private:

  }; //class OptEM
} // namespace Likelihood
