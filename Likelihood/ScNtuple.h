#ifndef ScNtuple_h
#define ScNtuple_h

#include "astro/SkyDir.h"

namespace Likelihood {

/** 
 * @class ScNtuple
 *
 * @brief NTuple class to represent spacecraft data.
 *
 * @author J. Chiang
 *    
 * $Header: */

class ScNtuple {
public:
   ScNtuple(){};
   ~ScNtuple(){};
   double time;
   astro::SkyDir zenDir;
   astro::SkyDir xAxis;
   astro::SkyDir zAxis;
   int inSaa;
};

} // namespace Likelihood
#endif // ScNtuple_h
