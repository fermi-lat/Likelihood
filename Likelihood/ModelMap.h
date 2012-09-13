/**
 * @file ModelMap.h
 * @brief Standalone class for computing a model map using a
 * BinnedLikelihood object.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/ModelMap.h,v 1.1 2012/09/13 19:24:29 jchiang Exp $
 */

#ifndef Likelihood_ModelMap_h
#define Likelihood_ModelMap_h

#include <string>
#include <vector>

namespace Likelihood {

class BinnedLikelihood;

/**
 * @class ModelMap
 * @brief Standalone class for computing a model map using a
 * BinnedLikelihood object.
 */

class ModelMap {

public:

   ModelMap(BinnedLikelihood & logLike);

   void writeOutputMap(const std::string & outfile,
                       std::string outtype="CMAP");

private:

   BinnedLikelihood & m_logLike;

   std::vector<float> m_outmap;

   void trimExtensions(const std::string & outfile,
                       const std::string & outtype);
};

} // namespace Likelihood

#endif // Likelihood_ModelMap_h
