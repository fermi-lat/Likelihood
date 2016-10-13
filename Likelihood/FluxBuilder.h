/**
 * @file FluxBuilder.h
 * @brief Builder class for creating flux-style xml files from 
 * Likelihood::Sources.
 * the flux package.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/FluxBuilder.h,v 1.4 2005/02/27 06:42:24 jchiang Exp $
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/FluxBuilder.h,v 1.4 2005/02/27 06:42:24 jchiang Exp $
 */

class FluxBuilder : public XmlBuilder {

public:

   FluxBuilder(double emin, double emax);

   virtual ~FluxBuilder();

   virtual void addSourceModel(SourceModel& srcModel);

   virtual void addSource(Source &src);
   
   virtual void write(std::string xmlFile);
      
private:

   typedef XERCES_CPP_NAMESPACE_QUALIFIER DOMElement DomElement;
   DomElement * m_srcLib;
   DomElement * m_allSrcsElt;
   void getSourceType(Source &src, std::string &srcType);
   DomElement * fluxSource(Source & src);
   DomElement * gammaSpectrum(optimizers::Function &spectrum);
   DomElement * srcDirection(optimizers::Function &dir);
   DomElement * solidAngle(double mincos, double maxcos);
   DomElement * galDiffuse(Source & src);
   DomElement * mapCubeSource(Source & src);

   std::vector<double> m_energies;
   void makeEnergyGrid(double emin, double emax, unsigned int nee=200);

   void addUnderscores(std::string &name);
};

} // namespace Likelihood

#endif // Likelihood_FluxBuilder_h
