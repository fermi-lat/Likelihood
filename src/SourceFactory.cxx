/** @file SourceFactory.cxx
 * @brief Implementation for the SourceFactory class, which applies the
 * Prototype pattern to return clones of various gamma-ray Sources.
 *
 * @author J. Chiang
 *
 * $Header$
 */

#include "Likelihood/SpectrumFactory.h"
#include "Likelihood/SourceFactory.h"

namespace Likelihood {

SourceFactory::SourceFactory() {
// Add a PointSource modeled by a PowerLaw as the default

// Note that the default constructor is used here, which means that
// exposure will not be computed.  A setDir(ra, dec, [true]) will
// cause the exposure to be computed and thus requires prior
// specification of the ROI cuts and spacecraft data.
   PointSource ptsrc;

// Add a nominal PowerLaw spectrum.  Note that one needs to reset the
// Parameters from the default and add sensible bounds.
   SpectrumFactory specFactory;
   Function *powerLaw = specFactory.makeFunction("PowerLaw");

// Use a nominal Parameter set for now with Prefactor = 10 (assuming a
// scaling of 1e-9, set below), Index = -2, and Scale = 100 (MeV).
// Set the bounds here as well. 
   std::vector<Parameter> params;
   powerLaw->getParams(params);
   params[0].setValue(10);            // Prefactor
   params[0].setScale(1e-9);
   params[0].setBounds(1e-3, 1e3);
   params[1].setValue(-2);            // Index
   params[1].setBounds(-3.5, -1);
   params[2].setValue(100);           // Scale (this is fixed by default)
   powerLaw->setParams(params);

   ptsrc.setSpectrum(powerLaw);

   addSource("PointSource", &ptsrc, true);
}

SourceFactory::~SourceFactory() {
   std::map<std::string, Source *>::iterator it = m_prototypes.begin();
   for (; it != m_prototypes.end(); it++)
      delete it->second;
}

void SourceFactory::addSource(const std::string &name, Source* src, 
                              bool fromClone) {
   if (!m_prototypes.count(name)) {
      if (fromClone) {
         m_prototypes[name] = src->clone();
      } else {
         m_prototypes[name] = src;
      }
   } else {
      std::cerr << "SourceFactory: A Source named "
                << name << " already exists!" << std::endl;
      assert(false);
   }
}

Source *SourceFactory::makeSource(const std::string &name) {
   assert(m_prototypes.count(name));
   return m_prototypes[name]->clone();
}

void SourceFactory::listSources() {
   std::cout << "SourceFactory Sources: " << std::endl;
   std::map<std::string, Source *>::const_iterator 
      it = m_prototypes.begin();
   for (; it != m_prototypes.end(); it++)
      std::cout << it->first << std::endl;
}

} // namespace Likelihood
