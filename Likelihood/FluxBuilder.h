/**
 * @file FluxBuilder.h
 * @brief Builder class for creating flux-style xml files from 
 * Likelihood::Sources.
 * the flux package.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/FluxBuilder.h,v 1.1 2004/02/18 21:13:25 jchiang Exp $
 */

#ifndef Likelihood_FluxBuilder_h
#define Likelihood_FluxBuilder_h

#include "Likelihood/XmlBuilder.h"

namespace optimizers {
   class Function;
}

namespace Likelihood {

/**
 * @class FluxBuilder
 * @brief This class provides methods for writing the source
 * information from Source objects as xml output appropriate for the
 * flux package.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/FluxBuilder.h,v 1.1 2004/02/18 21:13:25 jchiang Exp $
 */

class FluxBuilder : public XmlBuilder {

public:

   FluxBuilder();

   virtual ~FluxBuilder();

   virtual void addSource(Source &src);
   
   virtual void write(std::string xmlFile);
      
private:

   DomElement * m_srcLib;
   DomElement * m_allSrcsElt;

   void getSourceType(Source &src, std::string &srcType);
   DomElement * fluxSource(Source & src);
   DomElement * gammaSpectrum(optimizers::Function &spectrum);
   DomElement * srcDirection(optimizers::Function &dir);
   DomElement * solidAngle(double mincos, double maxcos);
   DomElement * galDiffuse(Source & src);

   std::vector<double> m_energies;
   void makeEnergyGrid(unsigned int nee=200);

   void addUnderscores(std::string &name);
};

} // namespace Likelihood

#endif // Likelihood_FluxBuilder_h
