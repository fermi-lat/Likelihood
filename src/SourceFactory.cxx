/** 
 * @file SourceFactory.cxx
 * @brief Implementation for the SourceFactory class, which applies the
 * Prototype pattern to return clones of various gamma-ray Sources.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SourceFactory.cxx,v 1.16 2003/07/18 16:21:23 jchiang Exp $
 */

#include <cassert>
#include <sstream>
#include "Likelihood/PointSource.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/SpatialMap.h"
#include "Likelihood/ConstantValue.h"
#include "Likelihood/SpectrumFactory.h"
#include "Likelihood/SourceFactory.h"
#include "Likelihood/LikelihoodException.h"

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

// Add the map-based Galactic Diffuse Emission model;
// assume that the FITS file is available in a standard place...
   std::string galfile;
   const char * root = ::getenv("LIKELIHOODROOT");
   if (!root) {  //use relative path from cmt directory
      galfile = "../src/test/Data/gas.cel";
   } else {
      galfile = std::string(root) + "/src/test/Data/gas.cel";
   }
   SpatialMap galacticModel(galfile);
   galacticModel.setParam("Prefactor", 1.1*pow(100., 1.1));

   try {
      DiffuseSource ourGalaxy(&galacticModel);
      ourGalaxy.setName("Milky Way");

// Provide ourGalaxy with a power-law spectrum.
      PowerLaw gal_pl(pow(100., -2.1), -2.1, 100.);
      gal_pl.setName("gal_pl");
      gal_pl.setParamScale("Prefactor", 1e-5);
      gal_pl.setParamTrueValue("Prefactor", pow(100., -2.1));
      gal_pl.setParamBounds("Prefactor", 1e-3, 1e3);
      gal_pl.setParamBounds("Index", -3.5, -1);

      ourGalaxy.setSpectrum(&gal_pl);

      addSource("Milky Way", &ourGalaxy, true);
   } catch (ParameterNotFound &eObj) {
      std::cerr << eObj.what() << std::endl;
      throw;
   } catch (LikelihoodException &likeException) {
      std::cerr << "Likelihood::SourceFactory: "
                << "Cannot create DiffuseSource Milkyway.\n"
                << likeException.what() << std::endl;
   }

// Add an extragalactic diffuse component.
   ConstantValue egNorm(1.);
   egNorm.setParam("Value", 1., false);   // fix to unity

   try {
      DiffuseSource extragalactic(&egNorm);
      extragalactic.setName("EG component");

      PowerLaw eg_pl(2.09e-3*pow(100., -2.1), -2.1, 100.);
      eg_pl.setName("eg_pl");
      eg_pl.setParamScale("Prefactor", 1e-7);
      eg_pl.setParamTrueValue("Prefactor", 2.09e-3*pow(100., -2.1));
      eg_pl.setParamBounds("Prefactor", 1e-5, 1e2);
      eg_pl.setParamBounds("Index", -3.5, -1);
      extragalactic.setSpectrum(&eg_pl);

      addSource("EG component", &extragalactic, true);
   } catch (ParameterNotFound &eObj) {
      std::cerr << eObj.what() << std::endl;
      throw;
   } catch (LikelihoodException &likeException) {
      std::cerr << "Likelihood::SourceFactory: "
                << "Cannot create DiffuseSource EG component.\n"
                << likeException.what() << std::endl;
   }
}

SourceFactory::~SourceFactory() {
   std::map<std::string, Source *>::iterator it = m_prototypes.begin();
   for (; it != m_prototypes.end(); it++)
      delete it->second;
}

void SourceFactory::addSource(const std::string &name, Source* src, 
                              bool fromClone) 
   throw(LikelihoodException) {
   if (!m_prototypes.count(name)) {
      if (fromClone) {
         m_prototypes[name] = src->clone();
      } else {
         m_prototypes[name] = src;
      }
   } else {
      std::ostringstream errorMessage;
      errorMessage << "SourceFactory::addSource: A Source named "
                   << name << " already exists.\n";
      throw LikelihoodException(errorMessage.str());
   }
}

void SourceFactory::replaceSource(Source* src, bool fromClone) {
   if (m_prototypes.count(src->getName())) {
      if (fromClone) {
         m_prototypes[src->getName()] = src->clone();
      } else {
         m_prototypes[src->getName()] = src;
      }
   } else {
      std::cerr << "SourceFactory::replaceSource: A Source named "
                << src->getName() << " does not yet exist.\n"
                << "Adding it instead. "
                << std::endl;
      addSource(src->getName(), src, fromClone);
   }
}

Source *SourceFactory::makeSource(const std::string &name) {
   assert(m_prototypes.count(name));
   return m_prototypes[name]->clone();
}

void SourceFactory::fetchSrcNames(std::vector<std::string> &srcNames) {
   if (!srcNames.empty()) srcNames.clear();
   std::map<std::string, Source *>::const_iterator it = m_prototypes.begin();
   for (; it != m_prototypes.end(); it++)
      srcNames.push_back(it->first);
}

} // namespace Likelihood
