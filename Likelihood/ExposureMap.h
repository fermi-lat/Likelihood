/**
 * @file ExposureMap.h
 * @brief ExposureMap class declaration.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef ExposureMap_h
#define ExposureMap_h

#include <vector>
#include <string>

namespace Likelihood {

/**
 * @class ExposureMap 
 *
 * @brief This is a Singleton class that encapulates and provides
 * exposure map information, primarily for use by the DiffuseSource
 * class for integrating the response functions over the spatial
 * distributions of those Sources.
 *
 * The exposure map can be read in from an existing file or computed
 * ab initio given the ROI cuts and spacecraft data.
 *
 * @author J. Chiang
 *
 * $Header$
 *
 */

class ExposureMap {

public:

   ~ExposureMap() {}

   void fetchRA(std::vector<double> &ra) {ra = s_ra;}
   void fetchDec(std::vector<double> &dec) {dec = s_dec;}
   void fetchExposure(std::vector<double> &exposure) {exposure = s_exposure;}

   static void readExposureFile(std::string exposureFile);

   static ExposureMap * instance();

protected:

   ExposureMap() {}

private:

   static ExposureMap * s_instance;

   static std::vector<double> s_ra;
   static std::vector<double> s_dec;
   static std::vector<double> s_exposure;
};

} // namespace Likelihood

#endif  // ExposureMap_h
