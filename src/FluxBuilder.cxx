/**
 * @file FluxBuilder.cxx
 * @brief Implementation for class to provide methods to write flux package
 * style xml files.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/FluxBuilder.cxx,v 1.6 2004/11/11 04:32:24 jchiang Exp $
 */

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>

#include "facilities/Util.h"

#include "optimizers/Dom.h"
#include "optimizers/Function.h"

#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/TrapQuad.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/SpatialMap.h"
#include "Likelihood/Source.h"
#include "Likelihood/FluxBuilder.h"

namespace Likelihood {

using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;
using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;

FluxBuilder::FluxBuilder() : XmlBuilder() {

   m_srcLib = optimizers::Dom::createElement(m_doc, "source_library");
   xml::Dom::addAttribute(m_srcLib, "title", "Likelihood_model");

   m_allSrcsElt = optimizers::Dom::createElement(m_doc, "source");

   makeEnergyGrid();
}

FluxBuilder::~FluxBuilder() {}

void FluxBuilder::addSource(Source & src) {
   DOMElement * srcElt = fluxSource(src);

   if (srcElt) {
      optimizers::Dom::appendChild(m_srcLib, srcElt);
      
      DOMElement * nestedSrc 
         = optimizers::Dom::createElement(m_doc, "nestedSource");
      std::string name = src.getName();
      addUnderscores(name);
      xml::Dom::addAttribute(nestedSrc, "sourceRef", name.c_str());
      optimizers::Dom::appendChild(m_allSrcsElt, nestedSrc);
   }
}

void FluxBuilder::write(std::string xmlFile) {
   facilities::Util::expandEnvVar(&xmlFile);

   std::string name = std::string("all_in_") 
      + facilities::Util::basename(xmlFile.c_str());
   xml::Dom::addAttribute(m_allSrcsElt, "name", name);
   optimizers::Dom::appendChild(m_srcLib, m_allSrcsElt);

   std::ofstream outFile(xmlFile.c_str());
   xml::Dom::prettyPrintElement(m_srcLib, outFile, std::string(""));
}

DOMElement * FluxBuilder::fluxSource(Source & src) {

   DOMElement * srcElt = optimizers::Dom::createElement(m_doc, "source");
   std::string name = src.getName();
   addUnderscores(name);
   xml::Dom::addAttribute(srcElt, "name", name);

   std::string sourceType;
   getSourceType(src, sourceType);

   Source::FuncMap & srcFuncs = src.getSrcFuncs();
   TrapQuad fluxIntegral(srcFuncs["Spectrum"]);
   if (sourceType == "PointSource" || sourceType == "Isotropic") {
      xml::Dom::addAttribute(srcElt, std::string("flux"),
                             fluxIntegral.integral(m_energies)/1e-4);
      DOMElement * specElt = gammaSpectrum(*srcFuncs["Spectrum"]);
      if (sourceType == "PointSource") {
         optimizers::Dom::appendChild(specElt, 
                                      srcDirection(*srcFuncs["Position"]));
      } else {
         optimizers::Dom::appendChild(specElt, solidAngle(-0.4, 1.0));
      }
      optimizers::Dom::appendChild(srcElt, specElt);
      return srcElt;
   } else if (sourceType == "GalDiffuse") {
      optimizers::Dom::appendChild(srcElt, galDiffuse(src));
      return srcElt;
   }
   return 0;
}

void FluxBuilder::getSourceType(Source &src, std::string & srcType) {
   Source::FuncMap & srcFuncs = src.getSrcFuncs();

   if (srcFuncs.count("Position")) {
      srcType = "PointSource";
   } else if (srcFuncs.count("SpatialDist") 
              && srcFuncs["SpatialDist"]->genericName() == "ConstantValue") {
      srcType = "Isotropic";
   } else if (srcFuncs.count("SpatialDist")
              && srcFuncs["SpatialDist"]->genericName() == "SpatialMap") {
      std::string fitsFile 
         = dynamic_cast<SpatialMap *>(srcFuncs["SpatialDist"])->fitsFile();
      std::string basename = facilities::Util::basename(fitsFile.c_str());
      if (basename == "gas.cel") {
         srcType = "GalDiffuse";
      }      
   } else {
      throw Exception(std::string("Likelihood::FluxBuilder::getSourceType:\n")
                      + "unknown source type");
   }
   return;
}   

DOMElement * FluxBuilder::gammaSpectrum(optimizers::Function & spectrum) {
   
   DOMElement * specElt = optimizers::Dom::createElement(m_doc, "spectrum");
   xml::Dom::addAttribute(specElt, std::string("escale"), std::string("MeV"));

   DOMElement * partElt = optimizers::Dom::createElement(m_doc, "particle");
   xml::Dom::addAttribute(partElt, std::string("name"), std::string("gamma"));

// Determine spectral type and set parameter values.
   DOMElement * spectralTypeElt 
      = optimizers::Dom::createElement(m_doc, "power_law");
   
   xml::Dom::addAttribute(spectralTypeElt, std::string("emin"), 
                          m_energies.front());
   xml::Dom::addAttribute(spectralTypeElt, std::string("emax"), 
                          m_energies.back());

// It might be better here to have the Function objects set their own
// parameter value attributes, but this keeps the coupling between
// Functions and the xml interface looser.
   if (spectrum.genericName() == "PowerLaw") {
      xml::Dom::addAttribute(spectralTypeElt, std::string("gamma"), 
                             -spectrum.getParamValue("Index"));
   } else if (spectrum.genericName() == "BrokenPowerLaw") {
      xml::Dom::addAttribute(spectralTypeElt, std::string("gamma"), 
                             -spectrum.getParamValue("Index1"));
      xml::Dom::addAttribute(spectralTypeElt, std::string("gamma2"), 
                             -spectrum.getParamValue("Index2"));
      xml::Dom::addAttribute(spectralTypeElt, std::string("ebreak"), 
                             spectrum.getParamValue("BreakValue"));
   }
   optimizers::Dom::appendChild(partElt, spectralTypeElt);
   optimizers::Dom::appendChild(specElt, partElt);
   return specElt;
}

DOMElement * FluxBuilder::srcDirection(optimizers::Function & dir) {
   DOMElement * dirElt = optimizers::Dom::createElement(m_doc, 
                                                        "celestial_dir");
   xml::Dom::addAttribute(dirElt, std::string("ra"), 
                          dynamic_cast<SkyDirFunction*>(&dir)->getDir().ra());
   xml::Dom::addAttribute(dirElt, std::string("dec"), 
                          dynamic_cast<SkyDirFunction*>(&dir)->getDir().dec());
   return dirElt;
}

DOMElement * FluxBuilder::solidAngle(double mincos, double maxcos) {
   DOMElement * solidAngle = optimizers::Dom::createElement(m_doc, 
                                                            "solid_angle");
   xml::Dom::addAttribute(solidAngle, std::string("mincos"), mincos);
   xml::Dom::addAttribute(solidAngle, std::string("maxcos"), maxcos);
   return solidAngle;
}

DOMElement * FluxBuilder::galDiffuse(Source & src) {
   Source::FuncMap & srcFuncs = src.getSrcFuncs();

   DOMElement * specElt = optimizers::Dom::createElement(m_doc, "spectrum");
   xml::Dom::addAttribute(specElt, std::string("escale"), std::string("MeV"));

   DOMElement * specClassElt = optimizers::Dom::createElement(m_doc,
                                                              "SpectrumClass");
   xml::Dom::addAttribute(specClassElt, std::string("name"), 
                          std::string("MapSpectrum"));

   if (srcFuncs["Spectrum"]->genericName() != "PowerLaw") {
      throw Exception(std::string("SourceModel::write_fluxXml:\n")
                      + "Galactic Diffuse spectral model is not a power-law.");
   } else {
// flux::MapSpectrum does not allow one to change the overall scaling
// of the FITS map used for modeling the Galactic diffuse emission. So
// we can only use the spectral index from our fit in this xml entry.
      std::ostringstream params;
      params << "100.," 
             << m_energies[0] << "," 
             << m_energies[m_energies.size()-1] << "," 
             << -srcFuncs["Spectrum"]->getParamValue("Index") << ","
             << "/sources/gas_gal.fits";
      xml::Dom::addAttribute(specClassElt, std::string("params"), 
                             params.str());
   }
   optimizers::Dom::appendChild(specElt, specClassElt);
   DOMElement * useSpecElt = optimizers::Dom::createElement(m_doc,
                                                            "use_spectrum");
   xml::Dom::addAttribute(useSpecElt, std::string("frame"), 
                          std::string("galaxy"));
   optimizers::Dom::appendChild(specElt, useSpecElt);
   return specElt;
}

void FluxBuilder::makeEnergyGrid(unsigned int nee) {
   RoiCuts * roiCuts = RoiCuts::instance();
   std::pair<double, double> elims = roiCuts->getEnergyCuts();
   double estep = log(elims.second/elims.first)/(nee-1);
   m_energies.reserve(nee);
   for (unsigned int i = 0; i < nee; i++) {
      m_energies.push_back(elims.first*exp(estep*i));
   }
}

void FluxBuilder::addUnderscores(std::string &name) {
// Replace spaces with underscores.
   std::replace(name.begin(), name.end(), ' ', '_');

// Prepend underscore if name starts with an integer character.
   if (static_cast<int>(*name.begin()) >= '0'
       && static_cast<int>(*name.begin()) <= '9') {
      name = "_" + name;
   }
}

} // namespace Likelihood
