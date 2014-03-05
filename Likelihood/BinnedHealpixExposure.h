/**
 * @file BinnedHealpixExposure.h
 * @brief All-Sky exposure map for use by SourceMap for DiffuseSource 
 * integrations
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/BinnedHealpixExposure.h,v 1.23 2014/02/19 20:30:35 jchiang Exp $
 */

#ifndef Likelihood_BinnedHealpixExposure_h
#define Likelihood_BinnedHealpixExposure_h

#include <stdexcept>
#include <string>
#include <vector>

#include "Likelihood/ExposureCube.h"
#include "Likelihood/BinnedExposure.h"

#include "evtbin/HealpixMap.h"

namespace st_app {
   class AppParGroup;
}

namespace astro {
   class SkyProj;
}

namespace Likelihood {

   class Observation;

/**
 * @class BinnedHealpixExposure
 * @brief This class encapsulates the calculation of and access to 
 * the integral of the effective area over live time.
 *
 * @author J. Chiang
 */

class BinnedHealpixExposure : public BinnedExposure {

public:

   BinnedHealpixExposure(const evtbin::HealpixMap & cmap, 
                  const Observation & observation, 
                  bool useEbounds=true,
                  const st_app::AppParGroup * pars=0);

   virtual ~BinnedHealpixExposure();

   virtual void writeOutput(const std::string & filename) const;

protected:

   void setMapGeometry(const evtbin::HealpixMap & cmap);

//private:

   const Observation * m_observation;

   std::vector<float> m_exposureMap;

   std::vector<double> m_energies;

   bool m_isGalactic;

   void computeHealpixMap(const evtbin::HealpixMap & cmap, const std::string filename);

};

} // namespace Likelihood

#endif // Likelihood_BinnedHealpixExposure_h
