/**
 * @file FluxBuilder.h
 * @brief Class to build flux-style XML files from Likelihood sources
 * @author J. Chiang
 */

#ifndef Likelihood_FluxBuilder_h
#define Likelihood_FluxBuilder_h

#include <string>

#include "xmlBase/rapidxml.hpp"

#include "Likelihood/XmlBuilder.h"

namespace Likelihood {

class Source;
class SourceModel;

class FluxBuilder : public XmlBuilder {
public:
   FluxBuilder(double emin, double emax);
   
   virtual ~FluxBuilder();
   
   void addSource(Source& src);
   
   void addSourceModel(SourceModel& srcModel);
   
   void write(std::string xmlFile);

private:
   double m_emin;
   double m_emax;
};

} // namespace Likelihood

#endif // Likelihood_FluxBuilder_h
