/**
 * @file FluxModel.h
 * @brief Class to encapsulate xml source models for the flux package.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef Likelihood_FluxModel_h
#define Likelihood_FluxModel_h

#include <string>

#include "xml/Dom.h"

namespace optimizers {
   class Function;
}

namespace Likelihood {

   class Source;

/**
 * @class FluxModel
 * @brief This class provides methods for writing the source
 * information from the Likelihood classes, especially SourceModel, as
 * xml output appropriate for the flux package.
 *
 * @author J. Chiang
 *
 * $Header$
 */

class FluxModel {

public:

   FluxModel();

   ~FluxModel();

   void addSource(Source &src);

   void write(std::string xmlFile);
      
private:

   xml::XmlParser * m_parser;
   DomDocument * m_doc;
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

#endif // Likelihood_FluxModel_h
