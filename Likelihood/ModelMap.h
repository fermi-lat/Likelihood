/**
 * @file ModelMap.h
 * @brief Standalone class for computing a model map using a
 * BinnedLikelihood object.
 * @author J. Chiang
 *
 * $Header$
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
                       const std::string & outtype="CMAP");

private:

   BinnedLikelihood & m_logLike;

   std::vector<float> m_outmap;

   void trimExtensions(const std::string & outfile,
                       const std::string & outtype);
};

} // namespace Likelihood

#endif // Likelihood_ModelMap_h
