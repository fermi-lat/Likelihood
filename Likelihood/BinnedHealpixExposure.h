/**
 * @file BinnedHealpixExposure.h
 * @brief All-Sky exposure map for use by SourceMap for DiffuseSource 
 * integrations
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/echarles/healpix_changes/Likelihood/Likelihood/BinnedHealpixExposure.h,v 1.4 2015/03/05 19:58:24 echarles Exp $
 */

#ifndef Likelihood_BinnedHealpixExposure_h
#define Likelihood_BinnedHealpixExposure_h

#include <stdexcept>
#include <string>
#include <vector>

#include "Likelihood/ExposureCube.h"
#include "Likelihood/BinnedExposureBase.h"

#include "evtbin/HealpixMap.h"
#include "healpix_map.h"

namespace st_app {
   class AppParGroup;
}

namespace astro {
  class HealpixProj;
}

namespace Likelihood {

   class Observation;

/**
 * @class BinnedHealpixExposure
 * @brief This class encapsulates the calculation of and access to 
 * the integral of the effective area over live time.
 *
 * @author E. Charles, J. Cohen-Tanugi
 */

class BinnedHealpixExposure : public BinnedExposureBase {

public:

   BinnedHealpixExposure(const evtbin::HealpixMap & cmap, 
			 const Observation & observation, 
			 bool useEbounds=true,
			 const st_app::AppParGroup * pars=0);

   BinnedHealpixExposure(const std::string& filename);   

   virtual ~BinnedHealpixExposure();

   virtual double operator()(double energy, double ra, double dec) const;

   virtual void writeOutput(const std::string & filename) const;

protected:

   void setMapGeometry(const evtbin::HealpixMap & cmap);

   void setMapGeometry(); 

//private:

   astro::HealpixProj* m_healpixProj;

   typedef Healpix_Map<float> ExposurePlane_t;
   std::vector<ExposurePlane_t> m_exposureMap;

   void computeHealpixMap(const evtbin::HealpixMap & cmap);

};

} // namespace Likelihood

#endif // Likelihood_BinnedHealpixExposure_h
