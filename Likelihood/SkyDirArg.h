/** 
 * @file SkyDirArg.h
 * @brief Declaration of SkyDirArg class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SkyDirArg.h,v 1.3 2003/07/19 04:38:02 jchiang Exp $
 */

#ifndef Likelihood_SkyDirArg_h
#define Likelihood_SkyDirArg_h

#include "optimizers/Arg.h"
#include "astro/SkyDir.h"

namespace Likelihood {

/** 
 * @class SkyDirArg
 *
 * @brief Concrete Arg subclass for encapsulating data of type astro::SkyDir
 *
 * @authors J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SkyDirArg.h,v 1.3 2003/07/19 04:38:02 jchiang Exp $
 */

class SkyDirArg : public optimizers::Arg {
    
public:
   
   SkyDirArg(astro::SkyDir dir) : m_val(dir) {}
   virtual ~SkyDirArg() {}

   void fetchValue(astro::SkyDir &dir) const {dir = m_val;}

private:

   astro::SkyDir m_val;

};

} // namespace Likelihood

#endif // Likelihood_SkyDirArg_h
